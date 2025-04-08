#include "Phages.h"

#include "Compute.h"
#include "Forces.h"
#include "grow.h"
#include "Integrate.h"
#include "Constants.h"
#include "UniformGrid.h"
#include "Cell.h"
#include "tools.h"

#include <omp.h>
#include <random>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

////////////////////////////////////////////////////////////////////////////////
// Math functions
////////////////////////////////////////////////////////////////////////////////

// This function returns random numbers ~U(min,max)
int RU_Numbers(int min, int max)
{
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(min,max); // This is the Uniform pdf
        int value = distrib(gen); // The Burst size is random now
        return value;
}

// This function returns random numbers ~N(mu,sigma)
double RN_Numbers(double mu, double sigma)
{
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<double> distrib{mu,sigma};; // This is the Normal pdf
        double value = distrib(gen); // The Burst size is random now
        return value;
}

// This function returns random numbers ~LogN(mu,sigma)
double RLogN_Numbers(double mu, double sigma)
{
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::lognormal_distribution<double> distrib{mu,sigma};; // This is the Log-Normal pdf
        double value = distrib(gen); // The Burst size is random now
        return value;
}

////////////////////////////////////////////////////////////////////////////////
// Phage specific functions
////////////////////////////////////////////////////////////////////////////////

// Function to define the source of phage
void SourcePhage(int targetID, Cell* cells, double dt)
{
        // Infect the cell and be sure to see it crash in short future
        cells[targetID].PhageCell = dt;
        cells[targetID].Shrink    = true;
        //std::cout << "Source cell: " << targetID << std::endl;
        
}

// Function to burst a cell
void BurstCell(int targetID, const int last_pos, Cell* cells, UniformGrid& Grid)
{
        UniformGrid::Address tgtAddress, topAddress;      // for storing addresses of cells in the Grid
        Cell& targetCell = cells[targetID];   // to be deleted
        Cell& topCell    = cells[last_pos];  // to be swapped
        
        // find address of targetCell and topCell
        tgtAddress = Grid.GetAddress(average(targetCell.Position));
        topAddress = Grid.GetAddress(average(topCell.Position));
        
        // send the topCell to the old targetID places
        cells[targetID] = topCell;
        Grid.Move(targetID,topAddress,tgtAddress);
        
        //std::cout << "Moving from top to the tgt" << std::endl;

        // send the targetCell to the top and remove it
        cells[last_pos] = targetCell;
        Grid.Move(last_pos,tgtAddress,topAddress);
        
        Grid.Remove(last_pos,topAddress);

        //std::cout << "Moving from tgt to the top" << std::endl;
        
}

// Function to update infection status (number & loads)
int UpdateInfection(int N_cells, Cell* cells, double dt)
{
        int infected = 0;
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                Cell& testCell = cells[cellID];
                if(testCell.PhageCell > 0)
                {
                        testCell.PhageCell = testCell.PhageCell + dt;
                        infected++;
                }
        }
        return infected;
}

// Function to obtain a list of cells that are affected
int getAffectedList(int cellID, Cell* cells, UniformGrid& Grid, int maxNeighbours, int** NeighbourList, int* infectList)
{
        Cell& testCell = cells[cellID];
        UniformGrid::Address Address = Grid.GetAddress(average(testCell.Position));
        int numNeighbours = Grid.GetNeighbours(cellID, Address, NeighbourList[cellID], maxNeighbours);
        numNeighbours = reduceNeighbours(testCell, cells, NeighbourList[cellID]);
        memcpy(infectList, NeighbourList[cellID], (numNeighbours)*sizeof(int));
        
        return numNeighbours;
}

// This function returns true/false if cell is close to phage
bool infectCell(DoubleCoord loc, int CellID, Cell* cells, double tolAds)
{
        DoubleCoord ct = getCellcenter(CellID, cells);
        double dist = sqrt(pow(ct.x-loc.x,2) + pow(ct.y-loc.y,2) + pow(ct.z-loc.z,2));

        if(dist < (tolAds) and cells[CellID].PhageCell == 0)
        {
                return true;
        }
        else
        {
                return false;
        }
}

////////////////////////////////////////////////////////////////////////////////
// Properties of Cells
////////////////////////////////////////////////////////////////////////////////

// This function returns the radial distance to (0,0,0) of any cell
double getCellRadius(int CellID, Cell* cells)
{
        DoubleCoord P       = cells[CellID].Position.p;
        DoubleCoord Q       = cells[CellID].Position.q;
        DoubleCoord center  = DoubleCoord(0.5*(Q.x+P.x), 0.5*(Q.y+P.y), 0.5*(Q.z+P.z));
        double radius       = sqrt(pow(center.x,2) + pow(center.y,2));
        
        return radius;
}

// This function returns the center coords. of any cell
double getCellangle(int CellID, Cell* cells)
{
        DoubleCoord center = getCellcenter(CellID,cells);
        
        double angle = atan(center.y / center.x);
        
        return angle;
}

// This function returns the center coords. of any cell
DoubleCoord getCellcenter(int CellID, Cell* cells)
{
        DoubleCoord P       = cells[CellID].Position.p;
        DoubleCoord Q       = cells[CellID].Position.q;
        DoubleCoord center  = DoubleCoord(0.5*(Q.x+P.x), 0.5*(Q.y+P.y), 0.5*(Q.z+P.z));
        
        return center;
}

// This function returns the number of cells inside the box where the loc coordinate belongs
int GetBoxNumber(UniformGrid& Grid, DoubleCoord loc)
{
        UniformGrid::Address Address;
        IntCoord id_loc;
        int num_cells;
        
        // Get the address of the loc
        Address = Grid.GetAddress(loc);
        
        // Get the indexes (x,y,z) of the Box that has the Address object
        id_loc = Grid.GetXY(Address);
        
        // Get the number of cells in the box
        num_cells = Grid.GetNumber(id_loc.x, id_loc.y, id_loc.z);
        
        return num_cells;
}

// This function returns a list of cells inside the box where the loc coordinate belongs
void GetBoxCell(UniformGrid& Grid, DoubleCoord loc, int* CellinBox)
{
        UniformGrid::Address Address;
        IntCoord id_loc;
        
        // Get the address of the loc
        Address = Grid.GetAddress(loc);
        
        // Get the indexes (x,y,z) of the Box that has the Address object
        id_loc = Grid.GetXY(Address);
        
        // Fill in the Cells in the box
        Grid.GetCellsinBox(id_loc.x, id_loc.y, id_loc.z, CellinBox);
}

// This function returns geometry features of the colony
double getColonyRadii(int N_cells, Cell* cells)
{
        double radii_val = 0.0;
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                double val = getCellRadius(cellID,cells);
                if(val >= radii_val)
                {
                        radii_val = val;
                }
        }
        return radii_val;
}

double getColonyHeight(int N_cells, Cell* cells)
{
        double height_val = 0.0;
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                DoubleCoord ct = getCellcenter(cellID, cells);
                if(ct.z >= height_val)
                {
                        height_val = ct.z;
                }
        }
        return height_val;
}

DoubleCoord getColonyCenter(int N_cells, Cell* cells)
{
        DoubleCoord centerloc(0.0,0.0,0.0);
        
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                DoubleCoord Ct = getCellcenter(cellID, cells);
                centerloc.x    = centerloc.x + Ct.x;
                centerloc.y    = centerloc.y + Ct.y;
                centerloc.z    = centerloc.z + Ct.z;
        }
        centerloc.x = centerloc.x / N_cells;
        centerloc.y = centerloc.y / N_cells;
        centerloc.z = centerloc.z / N_cells;
        return centerloc;
}

////////////////////////////////////////////////////////////////////////////////
// Here comes all the Cases PacMan, Border, Casper, etc
////////////////////////////////////////////////////////////////////////////////

// This function returns the IDs of cells at border
int getBorderCell(int N_cells, Cell* cells, int srcs, int* CellPhageList)
{
        double tol = 1.0*L_divide;
        
        double radii = getColonyRadii(N_cells, cells);
        
        std::cout << "Final Radii val: " << std::setprecision(9) << fixed << radii << std::endl;
        
        int aux = 0;
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                DoubleCoord center = getCellcenter(cellID,cells);
                double val = getCellRadius(cellID,cells);
                //if( abs(val - radii) < tol and aux < srcs and center.z < tol)
                if(val > 0.7*radii and aux < srcs and center.z < tol)
                {
                        CellPhageList[aux] = cellID;
                        cells[cellID].Resistant = false;
                        aux++;
                        std::cout << "Affected cell: " << cellID << " Radii val: " \
                                  << std::setprecision(9) << fixed << val \
                                  << " coords.: (" << center.x \
                                  << ","           << center.y \
                                  << ","           << center.z 
                                  << ")" << std::endl;
                }
        }
        return aux;
}

// This function returns the IDs of cells at center
int getCenterCell(int N_cells, Cell* cells, int srcs, int* CellPhageList)
{
        double tol = L_divide;
        
        double radii  = getColonyRadii(N_cells, cells);
        double height = getColonyHeight(N_cells, cells);
        
        DoubleCoord CenterColony = getColonyCenter(N_cells, cells);
        
        std::cout << "Colony center:: (" << CenterColony.x << " , " << CenterColony.y << " , " << CenterColony.z << ")" << std::endl;
        //std::cout << "Final Radii val: " << std::setprecision(9) << fixed << radii << std::endl;
        
        int aux = 0;
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                DoubleCoord center = getCellcenter(cellID,cells);
                double val = sqrt(pow(center.x-CenterColony.x,2.0) + pow(center.y-CenterColony.y,2.0));
                if( val < tol and aux < srcs and center.z > 0.6*height)
                {
                        CellPhageList[aux] = cellID;
                        aux++;
                        std::cout << "Affected cell: " << cellID << " Radii val: " \
                                  << std::setprecision(9) << fixed << val \
                                  << " coords.: (" << center.x \
                                  << ","           << center.y \
                                  << ","           << center.z 
                                  << ")" << std::endl;
                }
        }
        return aux;
}


// This function returns the IDs of cells at Plaque location (circle initial)
int getPlaqueCell(int N_cells, Cell* cells, int srcs, int* PlaqueList)
{
        double tol = L_divide;
        double min_r = 100.0;
        int index;
        
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                double val = getCellRadius(cellID,cells);
                if(val <= min_r)
                {
                        min_r = val;
                        index = cellID;
                }
        }
        
        DoubleCoord center = getCellcenter(index,cells);
        
        std::cout << "Smallest Radii val: " << std::setprecision(9) << fixed << min_r << std::endl;
        std::cout << "Affected cell: " << index << " Radii val: " \
                  << std::setprecision(9) << fixed << min_r \
                  << " coords.: (" << center.x \
                  << ","           << center.y \
                  << ","           << center.z 
                  << ")" << std::endl;
        
        PlaqueList[0] = index;
        return 1;
        
}

// This function returns the ID of cell at border, as PacMan
int getPacmanCell(int N_cells, Cell* cells, int srcs, int* CellPhageList)
{
        double tol = 2.0*L_divide;
        
        double radii = getColonyRadii(N_cells, cells);
        
        std::cout << "Final Radii val: " << std::setprecision(9) << fixed << radii << std::endl;
        
        int aux = 0;
        
        for(int cellID=0; cellID<N_cells; cellID++)
        {
                DoubleCoord center = getCellcenter(cellID,cells);
                double val   = getCellRadius(cellID,cells);
                if( abs(val - radii) < tol and aux < srcs and center.z < tol and center.y > tol and abs(center.x-center.y) < tol)
                {
                        CellPhageList[aux] = cellID;
                        aux++;
                        std::cout << "Affected cell: " << cellID << " Radii val: " \
                                  << std::setprecision(9) << fixed << val \
                                  << " coords.: (" << center.x \
                                  << ","           << center.y \
                                  << ","           << center.z 
                                  << ")" << std::endl;
                        break;
                }
        }
        return aux;
}



