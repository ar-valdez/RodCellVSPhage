#include "Simulation.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <random>
#include <fstream>
#include <vector>
#include <algorithm>    // std::unique, std::sort

#include <sys/stat.h>
#include <string>
#include <stdlib.h>

#include "Array.h"
#include "Cell.h"
#include "Cells.h"
#include "Phages.h"
#include "Compute.h"
#include "InputOutput.h"
#include "Constants.h"
#include "tools.h"
#include "UniformGrid.h"
#include "Nutrients.h"
#include "Neighbours.h"
#include "Forces.h"
#include "ClockIt.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

using namespace std;

//void SimSingleThread(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, \
//                     int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, \
//                     DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, \
//                     DoubleArray3D& Density2, DoubleArray2D& WallDensity, \
//                     DoubleArray2D& WallDensity1, DoubleArray2D& WallDensity2, \
//                     EnvArray3D& Environment, EnvArray3D& oldEnvironment, \
//                     AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, \
//                     AgaArray2D** oldFieldWall, CoordArray2D& Normal,const std::string DirName, \
//                     const int caseTest);


int GetProcessorCount()
{
#ifdef _WIN32
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);    
#endif
}

void RunSimulation(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, \
                   int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, \
                   DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, \
                   DoubleArray3D& Density2, DoubleArray2D& WallDensity, \
                   DoubleArray2D& WallDensity1, DoubleArray2D& WallDensity2, \
                   EnvArray3D& Environment, EnvArray3D& oldEnvironment, \
                   AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, \
                   AgaArray2D** oldFieldWall, CoordArray2D& Normal, const std::string DirName, \
                   const int caseTest)
{
        SimSingleThread(N_cells, old_cells, new_cells, NeighbourList, maxNeighbours, \
                        Grid, Files, append, Height, Density, Density1, Density2, \
                        WallDensity, WallDensity1, WallDensity2, Environment, oldEnvironment, \
                        FieldAgar, oldFieldAgar, FieldWall, oldFieldWall, Normal, DirName, \
                        caseTest);
}

////////////////////////////////////////////////////////////////////////////////
////////////////// Main function body //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void SimSingleThread(int N_cells, Cell* old_cells, Cell* new_cells, int** NeighbourList, \
                     int maxNeighbours, UniformGrid& Grid, OutputFiles Files, bool append, \
                     DoubleArray2D& Height, DoubleArray3D& Density, DoubleArray3D& Density1, \
                     DoubleArray3D& Density2, DoubleArray2D& WallDensity, \
                     DoubleArray2D& WallDensity1, DoubleArray2D& WallDensity2, \
                     EnvArray3D& Environment, EnvArray3D& oldEnvironment, \
                     AgaArray3D** FieldAgar, AgaArray3D** oldFieldAgar, AgaArray2D** FieldWall, \
                     AgaArray2D** oldFieldWall, CoordArray2D& Normal, const std::string DirName,
                     const int caseTest)
{
        // initialize time (generations)
        double dt = initial_dt;	// time step	
        double t = t0;			// current time
        int minx, maxx, miny, maxy, maxz;
        
        // File to write the number of cells and phage vs time
        stringstream filename;
        ofstream solFile;
        solFile.setf(ios::scientific, ios::floatfield);
        filename << DirName << "/Cell_time.dat";
        solFile.open(filename.str(), ios::out);
        
        // Files to store phage locations
        stringstream filePhage;
        ofstream outFile_Phage;
        //stringstream PhageFolder= DirName << "/Phage";
        //mkdir(PhageFolder.str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // Done only once
        
        int burst_count     = 0;
        int infected        = 0;
        int removed         = 0;
        int N_phages        = 0;
        double pressureCell = 0.0;

        // create place to store the props for phages
        std::vector<double> lx,ly,lz;
        std::vector<bool> mask;
        
	// **********************Initialize*******************************

        // counter to determine when to output data
        double NextOutTime                 = OutputTime;
        double NextUpdateTime              = UpdateTime;
        double HeightDensityNextUpdateTime = UpdateTime;
        double HeightDensityUpdateTime     = UpdateTime;
        
        bool OutFlag                 = true;
        bool UpdateFlag              = true;
        bool HeightDensityUpdateFlag = true;
        bool phage_case              = false; // Only for phage infection
        bool phage_source            = true; // only to seed phages
        bool debug                   = true;
        bool ThreeDcase              = false;

	// filter for doing some smoothing and averaging of fields
	DoubleArray2D Filter(FilterLen,FilterLen);
	double value = 1.0/(double)(FilterLen*FilterLen);

	for (int ii = 0; ii<FilterLen; ii++)
	{
		for (int jj = 0; jj<FilterLen; jj++)
		{
			Filter.Set(ii,jj,value);
		}
	}
	int FilterDim = FilterLen/2;

        DoubleArray3D Filter3D(FilterLen,FilterLen,FilterLen);
        double norm = 0;
        for (int ii=0; ii<FilterLen; ii++)
        {
                for (int jj=0; jj<FilterLen; jj++)
                {
                        for (int kk=0; kk<FilterLen; kk++)
                        {
                                value = exp(-(ii-FilterLen/2)*(ii-FilterLen/2)-(jj-FilterLen/2)*(jj-FilterLen/2)-(kk-FilterLen/2)*(kk-FilterLen/2));
                                norm = norm + value;
                        }
                }
        }
        for (int ii=0; ii<FilterLen; ii++)
        {
                for (int jj=0; jj<FilterLen; jj++)
                {
                        for (int kk=0; kk<FilterLen; kk++)
                        {
                                value = exp(-(ii-FilterLen/2)*(ii-FilterLen/2)-(jj-FilterLen/2)*(jj-FilterLen/2)-(kk-FilterLen/2)*(kk-FilterLen/2));
                                Filter3D.Set(ii,jj,kk,value/norm);
                        }
                }
        }

	// *******************Calculate important fields******************
        DoubleArray3D RoughDensity(Density.Size().x, Density.Size().y, Density.Size().z);
        DoubleArray3D RoughDensity1(Density1.Size().x, Density1.Size().y, Density1.Size().z);
        DoubleArray3D RoughDensity2(Density2.Size().x, Density2.Size().y, Density2.Size().z);
        DoubleArray3D insideColonyDen(Density.Size().x, Density.Size().y, Density.Size().z);

        DoubleArray2D RoughWallDensity(WallDensity.Size().x, WallDensity.Size().y);
        DoubleArray2D RoughWallDensity1(WallDensity1.Size().x, WallDensity1.Size().y);
        DoubleArray2D RoughWallDensity2(WallDensity2.Size().x, WallDensity2.Size().y);
        DoubleArray3D RoughDensityShiftP(Density.Size().x, Density.Size().y, Density.Size().z);
        DoubleArray3D RoughDensity1ShiftP(Density1.Size().x, Density1.Size().y, Density1.Size().z);
        DoubleArray3D RoughDensity2ShiftP(Density2.Size().x, Density2.Size().y, Density2.Size().z);
        DoubleArray2D RoughWallDensityShiftP(WallDensity.Size().x, WallDensity.Size().y);
        DoubleArray2D RoughWallDensity1ShiftP(WallDensity1.Size().x, WallDensity1.Size().y);
        DoubleArray2D RoughWallDensity2ShiftP(WallDensity2.Size().x, WallDensity2.Size().y);
        DoubleArray3D DensityShiftP(Density.Size().x, Density.Size().y, Density.Size().z);
        DoubleArray3D Density1ShiftP(Density1.Size().x, Density1.Size().y, Density1.Size().z);
        DoubleArray3D Density2ShiftP(Density2.Size().x, Density2.Size().y, Density2.Size().z);
        DoubleArray2D WallDensityShiftP(WallDensity.Size().x, WallDensity.Size().y);
        DoubleArray2D WallDensity1ShiftP(WallDensity1.Size().x, WallDensity1.Size().y);
        DoubleArray2D WallDensity2ShiftP(WallDensity2.Size().x, WallDensity2.Size().y);
        DoubleArray2D RoughHeight(Height.Size().x, Height.Size().y);
        
        RoughHeight.Initialize(cellRadius);
        RoughDensity.Initialize(0.0);
        RoughDensity1.Initialize(0.0);
        RoughDensity2.Initialize(0.0);
        insideColonyDen.Initialize(0.0);

        RoughWallDensity.Initialize(0.0);
        RoughWallDensity1.Initialize(0.0);
        RoughWallDensity2.Initialize(0.0);
        RoughDensityShiftP.Initialize(0.0);
        RoughDensity1ShiftP.Initialize(0.0);
        RoughDensity2ShiftP.Initialize(0.0);
        RoughWallDensityShiftP.Initialize(0.0);
        RoughWallDensity1ShiftP.Initialize(0.0);
        RoughWallDensity2ShiftP.Initialize(0.0);
        DensityShiftP.Initialize(0.0);
        Density1ShiftP.Initialize(0.0);
        Density2ShiftP.Initialize(0.0);
        WallDensityShiftP.Initialize(0.0);
        WallDensity1ShiftP.Initialize(0.0);
        WallDensity2ShiftP.Initialize(0.0);

	// Make rough wall
	DoubleArray2D Wall(BoxX,BoxY);
	for (int ii = 0; ii < BoxX; ii++)
	{
		for (int jj = 0; jj < BoxY; jj++)
		{
			Wall.Set(ii,jj,drand48()*wall_rough);
		}
	}
	std::cout << "Created wall" << std::endl;

	int* IDlist = new int[maxNeighbours];
	int IDlen = 0;

        // cell stress tensor
	Tensor stressTensor;

        // net force
	DoubleCoord Fnet;

        // dividing cells
	int* dividingCells = new int[maxCells];		// list of cells that need to divide at the end of a timestep
	
	// calculate fields, output to file
        CreateOutputFiles(0, Files, append);
        
        GetDensity(RoughDensity, RoughDensity1, RoughDensity2, insideColonyDen, \
                   Height, Grid, old_cells, minx, maxx, miny, maxy, maxz);
        ShiftDensity(RoughDensity, RoughDensity1, RoughDensity2, RoughWallDensity, \
                     RoughWallDensity1, RoughWallDensity2,RoughDensityShiftP, \
                     RoughDensity1ShiftP, RoughDensity2ShiftP,RoughWallDensityShiftP, \
                     RoughWallDensity1ShiftP, RoughWallDensity2ShiftP,BoxX, BoxY, BoxZ);
        RoughDensityShiftP.Output(Files.roughDensity,3);
        RoughDensity1ShiftP.Output(Files.roughDensity1,3);
        RoughDensity2ShiftP.Output(Files.roughDensity2,3);

        BoxAverage(RoughDensity, Density, RoughWallDensity, WallDensity, Filter3D, \
                   int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
        BoxAverage(RoughDensity1, Density1,RoughWallDensity1, WallDensity1, Filter3D, \
                   int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
        BoxAverage(RoughDensity2, Density2, RoughWallDensity2, WallDensity2,Filter3D, \
                   int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
        BoxAverage(RoughDensityShiftP, DensityShiftP, RoughWallDensityShiftP, WallDensityShiftP, \
                   Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
        BoxAverage(RoughDensity1ShiftP, Density1ShiftP,RoughWallDensity1ShiftP, WallDensity1ShiftP, \
                   Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
        BoxAverage(RoughDensity2ShiftP, Density2ShiftP, RoughWallDensity2ShiftP, WallDensity2ShiftP, \
                   Filter3D, int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
        DensityShiftP.Output(Files.density,3);
        Density1ShiftP.Output(Files.density1,3);
        Density2ShiftP.Output(Files.density2,3);
        WallDensityShiftP.Output(Files.walldensity);
        WallDensity1ShiftP.Output(Files.walldensity1);
        WallDensity2ShiftP.Output(Files.walldensity2);
        Height_Average(Height, RoughHeight, minx, maxx, miny, maxy);
        RoughHeight.Output(Files.roughheight);
        Smooth(RoughHeight, Height, Filter, FilterDim, minx, maxx, miny, maxy);
        Height.Output(Files.height);

        // find surface normal
        GetSurfaceNormal(Normal, Height, minx, maxx, miny, maxy);
        Normal.Output(Files.normal);

        int Nconv = 0;

        Nconv = UpdateEnvArray(&Environment, &oldEnvironment, FieldAgar, oldFieldAgar, \
                               FieldWall, oldFieldWall, DensityShiftP, Density1ShiftP, \
                               Density2ShiftP, WallDensityShiftP, WallDensity1ShiftP, \
                               WallDensity2ShiftP, minx, maxx, miny, maxy, maxz, Height, \
                               insideColonyDen);
    
	//ShiftGrowthrate(&Environment, &FieldWall,BoxX, BoxY, BoxZ);

	Environment.Output(Files.env,1);
        for (int level=0;level<maxLevels;level++)
        {
                FieldAgar[level]->Append(Files.aga,1);
                FieldWall[level]->Append(Files.wal);
        }
    
	// output cell data
	if (OutFlag)
	{
        	pressureCell = 0.0;
        	for (int cellID=0;cellID<N_cells;cellID++)
		{
			IntCoord XYAddress = Grid.GetXY( Grid.GetAddress(average(old_cells[cellID].Position)) );
			DoubleCoord T(0,0,0);
			F_surf_tension(old_cells[cellID], Grid, XYAddress, Height, Normal, Fnet, T);
			Output(Files.cells, cellID, t, old_cells[cellID], Fnet);
			
                        // update the pressure for the colony
                        double cell_press = getPressure(Fnet);
                        pressureCell = pressureCell + cell_press;
		}
	}
	
	fflush(Files.cells);
	OutFlag = false;
        CloseOutputFiles(Files);

	// **************************** Loop through time and evolve colony ****************************
        ClockIt     ti0, ti1, ti2, ti3, ti4, ti00, tiph;
        double      SNutr, Smovgro, Swrite, Sdivide, Supdateout, st, Srho, sph;
        SNutr=0; Smovgro=0; Swrite=0; Sdivide=0; Supdateout=0; st=0; Srho=0; sph=0;
        int koutput=0;

        rewind(Files.lineage);
        int max_simulSteps = (t_max - t) / dt;
        
        // Phage motion update coeffs
        double Dloc = Dph * 3600 * 1.0e+12;
        double Phageupdate = 2.0 * sqrt(Dloc * dt) / (3.0*L_divide*1.0e-04);
        
        if(ThreeDcase)
        {
                Phageupdate = 2.0 * sqrt(Dloc * dt) / (1.0e+03*3.0*L_divide*1.0e-04);
        }
        
        //std::cout << "Diff " << Dph << " Cte: " << Phageupdate << std::endl;
        
        // Phage adsorption coeffs
        double Tdiv   = 0.5 * pow(L_divide*1.0e-06,2.0) / (3600.0*Dph); //Time to move L_div in hours
        double tolAds = 0.75 * L_divide + L_divide / (60*Tdiv*Adsorption_rate);
        int jump      = 50;//min(int(latent_period / dt) , 50);
        
        //std::cout << "Tolerance ads.: " << tolAds << std::endl;
        //exit(0);
        
        //std::cout << t << " " << t0 << " " << initial_dt << " " << dt << std::endl;
        
        //exit(0);
       
        for(int iter = 0; iter < max_simulSteps; iter++)
        {

                // Compute the real time
                t = t + dt;

                if(N_cells > maxCells-1)
                {
                        std::cout << "Simulation ended. Max. number of cells reached" << std::endl;
                        break;
                }
                
		// update fields only when UpdateFlag is true
                ti0.start();
        
                if (UpdateFlag and ThreeDcase)
                {
                    UpdateFlag = false;
                    
                    // find density
                    GetDensity(RoughDensity, RoughDensity1, RoughDensity2, insideColonyDen, \
                               Height, Grid, old_cells, minx, maxx, miny, maxy, maxz);
                    ShiftDensity(RoughDensity, RoughDensity1, RoughDensity2, RoughWallDensity, \
                                 RoughWallDensity1, RoughWallDensity2,RoughDensityShiftP, \
                                 RoughDensity1ShiftP, RoughDensity2ShiftP,RoughWallDensityShiftP, \
                                 RoughWallDensity1ShiftP, RoughWallDensity2ShiftP,BoxX, BoxY, BoxZ);
                    BoxAverage(RoughDensity, Density, RoughWallDensity, WallDensity, Filter3D, \
                               int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
                    BoxAverage(RoughDensity1, Density1,RoughWallDensity1, WallDensity1, Filter3D, \
                               int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
                    BoxAverage(RoughDensity2, Density2, RoughWallDensity2, WallDensity2,Filter3D, \
                               int((Filter3D.m_Size.x)/2), minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
                    BoxAverage(RoughDensityShiftP, DensityShiftP, RoughWallDensityShiftP, \
                               WallDensityShiftP, Filter3D, int((Filter3D.m_Size.x)/2), \
                               minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
                    BoxAverage(RoughDensity1ShiftP, Density1ShiftP,RoughWallDensity1ShiftP, \
                               WallDensity1ShiftP, Filter3D, int((Filter3D.m_Size.x)/2), \
                               minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
                    BoxAverage(RoughDensity2ShiftP, Density2ShiftP, RoughWallDensity2ShiftP, \
                               WallDensity2ShiftP,Filter3D, int((Filter3D.m_Size.x)/2), \
                               minx, maxx, miny, maxy, maxz,BoxX,BoxY,BoxZ);
                    Height_Average(Height, RoughHeight, minx, maxx, miny, maxy);
                    //			RoughHeight.Output(Files.height);
                    Smooth(RoughHeight, Height, Filter, FilterDim, minx, maxx, miny, maxy);
                    //			Height.Output(Files.height);
                    
                    // find surface normal
                    GetSurfaceNormal(Normal, Height, minx, maxx, miny, maxy);
                    
                   // find nutrient concentrations
                    Nconv = UpdateEnvArray(&Environment, &oldEnvironment, FieldAgar, oldFieldAgar, \
                                           FieldWall, oldFieldWall, DensityShiftP, Density1ShiftP, \
                                           Density2ShiftP, WallDensityShiftP, WallDensity1ShiftP, \
                                           WallDensity2ShiftP, minx, maxx, miny, maxy, maxz, Height, \
                                           insideColonyDen);
                    // ShiftGrowthrate(&Environment, &FieldWall,BoxX, BoxY, BoxZ);
                    if (Nconv>=maxIter)
                    {
                        std::cout << "Environment has not converged" << std::endl;
                    }
                }
                SNutr+=ti0.stop()/1000.0;
                ti00.start();
                if (HeightDensityUpdateFlag and ThreeDcase)
                {
                    HeightDensityUpdateFlag = false;
	            
                    // find height
                    GetHeight(RoughDensity, Height, Grid, old_cells, minx, maxx, miny, maxy, maxz);
                    Height_Average(Height, RoughHeight, minx, maxx, miny, maxy);
                    Smooth(RoughHeight, Height, Filter, FilterDim, minx, maxx, miny, maxy);
                    // find surface normal
                    GetSurfaceNormal(Normal, Height, minx, maxx, miny, maxy);
                }
                Srho+=ti00.stop()/1000.0;
                
                ti1.start();
                if (OutFlag)
                {


                        CreateOutputFiles(koutput++, Files, append);
                        pressureCell = 0.0;
                        for (int cellID=0;cellID<N_cells;cellID++)
                        {
                                // output cell data
                                mean_stress(old_cells[cellID], old_cells, NeighbourList[cellID], \
                                            caseTest,Grid, Wall, Height, Normal, stressTensor, Fnet);
                                Output(Files.cells, cellID, t, old_cells[cellID], stressTensor);

                                // update the pressure for the colony
                                double cell_press = getPressure(Fnet);
                                pressureCell = pressureCell + cell_press;
                        }
                }

                // For debug
                //std::cout << "Before move and grow" << std::endl << std::endl;

                // Update neighbour lists
                getNeighbours(old_cells, N_cells, Grid, NeighbourList, maxNeighbours);

                // calculate forces and move cells
                int numDivide = 0; // number of cells dividing
                int cellID;
                #pragma omp parallel for default(shared) private(cellID) schedule(static)
                for(cellID=0;cellID<N_cells;cellID++)
                {

                        GrowCell(old_cells[cellID], cellID, dt, dividingCells, numDivide, Environment, FieldWall, Grid, ThreeDcase);
                        //std::cout << "Growth done good" << std::endl;
                        
                        MoveCell(cellID, Grid, old_cells, new_cells, NeighbourList[cellID], dt, Height, Normal, Wall, caseTest);
                        //std::cout << "Move done good" << std::endl;
                }
	        //	fflush(Files.cells);
                
                // For debug
                //std::cout << "After move and grow" << std::endl << std::endl;

                Smovgro+=ti1.stop()/1000.0;

                ti3.start();

                // For debug
                //std::cout << "Before divide" << std::endl << std::endl;
	        
	        // Division (must be done after integration step because new neighbours are created)
	        for(int cellCount = 0; cellCount < numDivide; cellCount++)
	        {
		        // ID of cell that is undergoing division
		        const int cellID = dividingCells[cellCount];

		        DivideCell(cellID, N_cells, new_cells, Grid, NeighbourList[cellID], Wall, Height, Normal, t);
		        fprintf(Files.lineage, "%d %d\n", cellID, N_cells);
		        N_cells++;

		        // make list of all neighbours, new cell, and old cell
		        IDlen = NeighbourList[cellID][0];
		        memcpy(IDlist, NeighbourList[cellID], (IDlen+1)*sizeof(int)); // copy neighbours
		        IDlist[0] = cellID;
		        IDlist[IDlen+1] = N_cells-1;
		        IDlen+=2;

		        // update the neighbours of all of these cells
		        getNeighbours(new_cells, Grid, NeighbourList, maxNeighbours, IDlist, IDlen);

                        if(N_cells > maxCells-1)
                        {
                                std::cout << "Simulation ended. Max. number of cells reached" << std::endl;
                                break;
                        }
	        }

                // For debug
                //std::cout << "After divide" << std::endl << std::endl;
                Sdivide+=ti3.stop()/1000.0;

////////////////////////////////////////////////////////////////////////////////
/////////////////////// Traditional tests //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
                if(caseTest == 1 and phage_source)
                {
                        std::cout << "Make a phage-source-free case" << std::endl;
                        phage_source = false;
                        phage_case   = false;
                }
////////////////////////////////////////////////////////////////////////////////
//////////////// START Phage dynamics equations ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
                // Start the phage infection
                tiph.start();
                
                if(phage_source and t >= colony_time)
                {
                        std::cout << std::endl << "Inserting the sources of phage" << std::endl;
                        phage_source = false;
                        phage_case   = true;
                        int InfectSize = int(Infect_frac * colony_size / 100.0);
                        
                        // Prepare and fill in the lists of cells to infect
                        int* toPhageList = new int[InfectSize];
                        
                        if(caseTest == 2)
                        {
                                std::cout << "Make the PacMan Case" << std::endl;
                                InfectSize = getPacmanCell(N_cells, new_cells, InfectSize, toPhageList);
                        }
                        else if(caseTest == 3)
                        {
                                std::cout << "Make the Cells in border" << std::endl;
                                InfectSize = getBorderCell(N_cells, new_cells, InfectSize, toPhageList);
                        }
                        else if(caseTest == 11)
                        {
                                std::cout << "Make infection random on top Case" << std::endl;
                                for(int k=0; k < InfectSize; k++)
                                {
                                        SourcePhage(N_cells-k-1, new_cells,dt);
                                }
                                InfectSize = 0;
                        }
                        else if(caseTest == 12)
                        {
                                std::cout << "Inserting the sources of phage by Uniform placement" << std::endl;
                        
                                
                                double* pos_x = new double[InfectSize];
                                double* pos_y = new double[InfectSize];
                                double* pos_z = new double[InfectSize];
                                double radii  = getColonyRadii(N_cells, new_cells);
                                double height = getColonyHeight(N_cells, new_cells);
                                
                                
                                for(int k=0; k<InfectSize; k++)
                                {
                                        double PI = 4.0*atan(1.0);
                                        double r = drand48() * 0.01 * radii;
                                        double h = (1.0 - 0.01*drand48()) * height;
                                        double theta = (2.0*PI)*drand48();
                                        double psi = (1.0*PI)*drand48();
                                        double cx = r*sin(psi)*cos(theta);
                                        double cy = r*sin(psi)*sin(theta);
                                        double cz = h*cos(psi);
                                        lx.push_back(cx);
                                        ly.push_back(cy);
                                        lz.push_back(cz);
                                        
                                        mask.push_back(true);
                                }
                                N_phages = mask.size();
                                
                                delete[] pos_x;
                                delete[] pos_y;
                                delete[] pos_z;
                                InfectSize = 0;
                        }
                        else if(caseTest == 13)
                        {
                                std::cout << "Make the Cells in Center" << std::endl;
                                InfectSize = getCenterCell(N_cells, new_cells, InfectSize, toPhageList);
                        }
                        else
                        {
                                std::cout << "Phage case not implemented choose [2-13]" << std::endl;
                                std::cout << "caseTest :: 2 - PacMan" << std::endl;
                                std::cout << "caseTest :: 3 - Border Cells" << std::endl;
                                std::cout << "caseTest :: 11 - Cells random top" << std::endl;
                                std::cout << "caseTest :: 12 - Phage random uniform" << std::endl;
                                std::cout << "caseTest :: 13 - Center Cells" << std::endl;
                                InfectSize = 0;
                                exit(0);
                        }
                        
                        for(int k=0; k < InfectSize; k++)
                        {
                                int loc = toPhageList[k];
                                SourcePhage(loc, new_cells,dt);
                        }
                        delete[] toPhageList;
                        
                        std::cout << "Done Inserting the sources of phage" << std::endl << std::endl;
                }

                // Exit the simulation if extinction
                if(N_cells <= 2 and koutput > 10 and phage_source == false)
                {
                        std::cout << "Simulation ended. Colony extincted" << std::endl;
                        std::cout << N_cells << std::endl;
                        break;
                }
                
                // For debug
                //std::cout << "Before phage part" << std::endl << std::endl;
                
                if(phage_source == false and phage_case and iter % jump == 0)
                {
                        // Get the number of infected cells and update their loads
                        infected = UpdateInfection(N_cells, new_cells, dt * jump);
                        
                        // Get the IDs of infected cells ready to burst
                        std::vector<int> deadList;
                        
                        for(int k=0; k<N_cells; k++)
                        {
                                if(new_cells[k].PhageCell >= latent_period and new_cells[k].PhageCell < 10.0)
                                {
                                        deadList.push_back(k);
                                }
                        }

                        //Sort cells in the deadList and remove duplicates
                        std::sort(deadList.begin(), deadList.end());
                        deadList.erase(std::unique(deadList.begin(), deadList.end()), deadList.end());
                        
                        // Proceed to expand phage population and remove bursted cells
                        for(int infec_cellID=0; infec_cellID<deadList.size(); infec_cellID++)
                        {
                                // Get index on cell list
                                int loc_tgt = deadList[infec_cellID];
                                
                                // Store the location of cell to burst
                                DoubleCoord ct = getCellcenter(loc_tgt, new_cells);
                                
                                // Get Abortive status of cell
                                bool abortive = new_cells[loc_tgt].Abortive;

                                // Get Resistant status of cell
                                bool resistant = new_cells[loc_tgt].Resistant;
                                
                                // Knowing the location of bursted cell release random phage there
                                if(resistant == false and abortive==false and new_cells[loc_tgt].PhageCell < 10.0)
                                {
                                        for(int k=0; k<ProlePhage; k++)
                                        {
                                                // Place phages in new random places
                                                double PI    = 4.0*atan(1.0);
                                                double r     = drand48()*0.1;
                                                double theta = (2.0*PI)*drand48();
                                                double psi   = (1.0*PI)*drand48();

                                                lx.push_back(ct.x + r*sin(psi)*cos(theta));
                                                ly.push_back(ct.y + r*sin(psi)*sin(theta));
                                                lz.push_back(ct.z + 0.0*r*cos(psi));

                                                mask.push_back(true);
                                        }
                                        N_phages = N_phages + ProlePhage;
                                }
                                if(resistant == false)
                                {
                                        // Remove the loc_tgt cell from all the Neighbors
                                        for(int k=1; k< NeighbourList[loc_tgt][0]; k++)
                                        {
                                                int pos = NeighbourList[loc_tgt][k];
                                                removeKeyNeighbours(NeighbourList[pos], NeighbourList[pos][0], loc_tgt);
                                        }
                                        // Now we proceed to burst the cell
                                        BurstCell(loc_tgt, N_cells-1, new_cells, Grid);
                                        burst_count++;
                                        N_cells--;
                                }
                                
                                // For debug
                                //std::cout << "Killed stuff ID:: " << infec_cellID + 1 << std::endl;
                        }
                        deadList.clear();
                        
                        // Update NeighbourList array
                        //getNeighbours(new_cells, N_cells, Grid, NeighbourList, maxNeighbours);

                        //Now we infect by proximity to phage location
                        int ph;
                        #pragma omp parallel for default(shared) private(ph) schedule(static)
                        for(ph=0; ph<mask.size(); ph++)
                        {
                                #pragma omp critical
                                {
                                if(mask[ph])
                                        {
                                                DoubleCoord pos = DoubleCoord(lx[ph],ly[ph],lz[ph]);
                                                int count = GetBoxNumber(Grid, pos);
                                                if(count > 1)
                                                {
                                                        int* CellinBox = new int[count];
                                                        GetBoxCell(Grid, pos, CellinBox);
                                                        for(int k=0; k<count; k++)
                                                        {
                                                                int IDcell = CellinBox[k];
                                                                bool infect = infectCell(pos, IDcell, new_cells, tolAds);
                                                                bool resistant = new_cells[IDcell].Resistant;
                                                                if(infect and resistant==false)
                                                                {
                                                                        new_cells[IDcell].PhageCell = dt;
                                                                        new_cells[IDcell].Shrink    = true;
                                                                        mask[ph] = false;
                                                                        removed++;
                                                                        break;
                                                                }
                                                                else if(infect and resistant==true)
                                                                {
                                                                        new_cells[IDcell].PhageCell = 0.0;
                                                                        new_cells[IDcell].Shrink    = true;
                                                                }
                                                        }
                                                delete[] CellinBox;
                                                }
                                                // Now we move only the free phage
                                                if(mask[ph])
                                                {
                                                        double PI = 4.0*atan(1.0);
                                                        double r = drand48();
                                                        double theta = (2.0*PI)*drand48();
                                                        double psi = (1.0*PI)*drand48();
                                                        lx[ph] = lx[ph] + r*sin(psi)*cos(theta) * Phageupdate * jump;
                                                        ly[ph] = ly[ph] + r*sin(psi)*sin(theta) * Phageupdate * jump;
                                                        lz[ph] = lz[ph]; //+ r*cos(psi) * Phageupdate * (BoxZ+BoxZAgar) * jump / (BoxX);
                                                }
                                        }
                                }
                        }
                }

                // For debug
                //std::cout << "After phage part" << std::endl << std::endl;
                sph+=tiph.stop()/1000;
                                
////////////////////////////////////////////////////////////////////////////////
//////////////// END Phage dynamics equations //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

                ti2.start();
                
	        // output fields to file
	        if (OutFlag)
	        {
		        OutFlag = false;

		        // save restart file with cell information
		        SaveCells(Files.restart, new_cells, N_cells, t, dt);

                        if(mask.size() > 0)
                        {
                                outFile_Phage.setf(ios::scientific, ios::floatfield);
                                
                                filePhage.str("");
                                filePhage << DirName << "/Phage/phage_" << int(koutput) << ".dat";
                                
                                outFile_Phage.open(filePhage.str(), ios::out);
                                for(int k=0; k<mask.size(); k++)
                                {
                                        if(mask[k])
                                        {
                                                outFile_Phage << std::setprecision(4) << t << " " \
                                                              << k     << " " \
                                                              << lx[k] << " " \
                                                              << ly[k] << " " \
                                                              << lz[k] << std::endl;
                                                
                                        }
                                }
                                outFile_Phage.close();
                        }

                        if(debug)
                        {
                                std::cout << std::setprecision(2) << fixed << "t = " << 100.0 * t/t_max << " %;\t" << N_cells\
                                          << " living cells, " << infected << " infected cells, " \
                                          << burst_count << " killed cells, " << N_phages - removed << " phages. Sim. Time = " \
                                          << int(st/60.0/60.0) << ":" << int(st/60.0)%60 << ":" << int(st)%60 \
                                          << "; File written: " << koutput \
                                          << ". Time cost = " << SNutr << " - " << Srho << " - " << Smovgro \
                                          << " - " << Swrite << " - " << Sdivide << " - " << Supdateout << " - " << sph << std::endl;
                        }
                        else
                        {
                                std::cout << std::setprecision(2) << fixed << "t = " << 100.0 * t/t_max << " %;\t" << N_cells \
                                          << " living cells, " << infected << " infected cells, " \
                                          << burst_count << " killed cells, " << N_phages - removed  << " phages. Sim. Time = " \
                                          << int(st/60.0/60.0) << ":" << int(st/60.0)%60 << ":" << int(st)%60 \
                                          << "; File written: " << koutput << std::endl;
                        }

                        // Update NeighbourList array
                        //getNeighbours(new_cells, N_cells, Grid, NeighbourList, maxNeighbours);
                        
                        // Evaluate geometry features
                        double radii_val = getColonyRadii(N_cells, new_cells);
                        double height_val = getColonyHeight(N_cells, new_cells);
                        
                        solFile << std::setprecision(4) << t << " " \
                                << N_cells     << " " << burst_count << " "\
                                << infected << " " << pressureCell << " "\
                                << N_phages - removed  << " "\
                                << radii_val << " " << height_val << std::endl;

                        st += SNutr+Smovgro+Swrite+Sdivide+Supdateout+Srho+sph;
                        SNutr=0; Smovgro=0; Swrite=0; Sdivide=0; Supdateout=0; Srho=0; sph=0;
                        fflush(stdout);
                        CloseOutputFiles(Files);
                        
	        }
                
                Swrite+=ti2.stop()/1000.0;
                
                ti4.start();

	        // switch positions of old and new cells
	        {
	                Cell* temp = new_cells;
	                new_cells  = old_cells;
	                old_cells  = temp;
	        }
	        
////////////////////////////////////////////////////////////////////////////////
////////////////// determine if we're writing output next time /////////////////
////////////////////////////////////////////////////////////////////////////////

	        OutFlag = (NextOutTime<=0);
	        NextOutTime = (OutFlag ? OutputTime : NextOutTime) - dt;

	        // determine if we're updating fields
	        UpdateFlag = (NextUpdateTime<=0);
	        NextUpdateTime = (UpdateFlag ? UpdateTime: NextUpdateTime) - dt;
                
                HeightDensityUpdateFlag = (HeightDensityNextUpdateTime<=0);
                HeightDensityNextUpdateTime = (HeightDensityUpdateFlag? HeightDensityUpdateTime: HeightDensityNextUpdateTime) - dt;

                Supdateout+=ti4.stop()/1000.0;

        }

        fflush(Files.lineage);
        
        solFile.close();
        
        lx.clear();
        ly.clear();
        lz.clear();
        mask.clear();
        std::cout << "Done" << std::endl;
}

