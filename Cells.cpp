#include "Cells.h"

#include "Compute.h"
#include "Forces.h"
#include "grow.h"
#include "Integrate.h"
#include "Constants.h"
#include "UniformGrid.h"
#include "Cell.h"
#include "Phages.h"

#include <omp.h>
#include <random>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

// This function returns the pressure in a single cell
double getPressure(const Tensor T)
{
        double p = (T.xx + T.yy + T.zz) / 3.0;
        return p;
}

// This function returns the pressure in a single cell
double getPressure(const DoubleCoord F)
{
        double p = (F.x + F.y + F.z) / 3.0;
        return p;
}

// Main function to move cell
void MoveCell(int cellID, UniformGrid& Grid, const Cell* old_cells, Cell* new_cells, \
              const int* neighbours, double dt, DoubleArray2D& Height, \
              CoordArray2D& Normal, DoubleArray2D& Wall, const int caseTest)
{
        UniformGrid::Address oldAddress, newAddress;	// for storing addresses of cells in the Grid

        DoubleCoord v, Fnet;
        IntCoord XYAddress;
        
        const Cell& oldCell = old_cells[cellID];
        Cell& newCell       = new_cells[cellID];

        // gives the current neighbours of the cell
        oldAddress = Grid.GetAddress(average(oldCell.Position));
        XYAddress  = Grid.GetXY(oldAddress);
        
        //std::cout << "Ready Move cell ID:" << cellID << std::endl;

        // integrates one step, updates positions from old to new
        integrate(dt, cellID, old_cells, new_cells, neighbours, Height, Normal, Grid, XYAddress, Wall, caseTest);

        //std::cout << "Integrated cell ID:" << cellID << std::endl;

        // check if the cell has moved out of its box
        newAddress = Grid.GetAddress(average(newCell.Position));
        if (newAddress.a!=oldAddress.a)
        {
                #pragma omp critical
                {
	                //std::cout << "Move? cell ID:" << cellID << std::endl;
	                Grid.Move(cellID, oldAddress, newAddress);
	                //std::cout << "Moved cell ID:" << cellID << std::endl;
                }
        }
        //std::cout << "Move done good" << std::endl;
        
        // Use the number of Neighbors to modulate growthrate
        new_cells[cellID].Nb = neighbours[0] + 1;
}

bool CheckDivide(double DivCrit, double Lcurrent, double Loc_Ldiv)
{
        double val = RN_Numbers(3.5, 1);
        double sfbet = val > 3 ? val : 3;
        double P = 0.5 + (1.0/PI) * atan( sfbet * (Lcurrent - Loc_Ldiv));
        
        if (P >= DivCrit)
        {
                return true;
        }
        else
        {
                return false;
        }
}

// Main function to grow cell
void GrowCell(Cell& cell, int cellID, double dt, int* dividingCells, int& numDivide, \
              EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid, bool Nutrients)
{
        if(cell.PhageCell == 0 and cell.Length < 1.05*cell.Ldiv)
        {
                grow(dt, cell, Env, Wal, Grid, Nutrients);
                int index;
                double DivCriteria = RN_Numbers(0.6,0.1);
                // check if cell will divide
                bool div_status = CheckDivide(DivCriteria, cell.Length, cell.Ldiv);
                //if (cell.Length > L_divide) // Old division criteria
                
                if(cell.Length < 0.6*cell.Ldiv)
                {
                        div_status = false;
                }
                
                if(cell.Shrink)
                {
                        div_status = false;
                }
                
                if(abs(cell.Position.q.z - cell.Position.p.z) > 1.5*cell.Radius )
                {
                        div_status = false;
                }
                
                
                if(div_status)
                {
                        #pragma omp critical
                        {
                                index = numDivide++;
                        }
                         dividingCells[index] = cellID;
                }
                //std::cout << "Growth done good" << std::endl;
        }
        
        // Update the age of the cell
        cell.AgeCell = cell.AgeCell + dt;
}

// Main function to divide cell
void DivideCell(int parentID, int daughterID, Cell* cells, UniformGrid& Grid, \
                const int* neighbours, DoubleArray2D& Wall, DoubleArray2D& Height, \
                CoordArray2D& Normal, double t)
{
	if(cells[parentID].PhageCell == 0)
	{
	        UniformGrid::Address oldAddress, newAddress;	// for storing addresses of cells in the Grid
	        //double F_centre;
	        DoubleCoord Fnet;
	        Tensor stressTensor;
	        Cell& parentCell = cells[parentID];
	        Cell& daughterCell = cells[daughterID];

	        // find address on the mother cell
	        oldAddress = Grid.GetAddress(average(parentCell.Position));
	        
	        // remove the ID from the grid
	        Grid.Remove(parentID, oldAddress);

	        // divide and create a new cell with ID N_cells
	        divide(parentCell, daughterCell, t);

	        // add parent and daughter to grid
	        Grid.Add(parentID, Grid.GetAddress(average(parentCell.Position)));
	        Grid.Add(daughterID, Grid.GetAddress(average(daughterCell.Position)));
	}
}

