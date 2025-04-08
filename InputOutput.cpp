#include "InputOutput.h" 

#include <sys/stat.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdlib.h>

#include <chrono>
#include <random>

#include "Cell.h"
#include "Constants.h"
#include "Forces.h"
#include "Phages.h"

using namespace std;

// input and output functions
// also initializes the starting positions of cells
void CreateOutputFileLineage(int OutputID, OutputFiles& Files, bool append)
{
        // create output file lineage
        char lineage_name[500];

        // concatenate filenames with suffix
        strcpy(lineage_name,DirName);
        strcat(lineage_name,"/lineage");
        mkdir(lineage_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(lineage_name,"%s/%d",lineage_name,OutputID);
        strcat(lineage_name,".dat");

        // open files for output

        Files.lineage = fopen(lineage_name, "w");	// file to store lineage
        if (Files.lineage == NULL) {
                std::cout << "Can't open lineage file." << std::endl;
                exit(1);
        }
}

void CloseOutputFileLineage(OutputFiles& Files)
{
    fclose(Files.lineage);
}

void CreateOutputFiles(int OutputID, OutputFiles& Files, bool append)
{
	// create output files
        char cell_name[500], restart_name[500], roughDensity_name[500], roughDensity1_name[500];
        char roughDensity2_name[500], density_name[500], density1_name[500], density2_name[500];
        char walldensity_name[500], walldensity1_name[500], walldensity2_name[500];
        char roughHeight_name[500], height_name[500], normal_name[500], env_name[500];
        char aga_name[500], wal_name[500];
        char phage_name[500];

	// concatenate filenames with suffix
        strcpy(cell_name,DirName);
        strcat(cell_name,"/Cells");
        mkdir(cell_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(cell_name,"%s/%d",cell_name,OutputID);
        strcat(cell_name,".dat");

        strcpy(phage_name,DirName);
        strcat(phage_name,"/Phage");
        mkdir(phage_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        //sprintf(phage_name,"%s/%d",phage_name,OutputID);
        //strcat(phage_name,".dat");
        
        strcpy(restart_name,DirName);
        strcat(restart_name,"/Restart");
        mkdir(restart_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(restart_name,"%s/%d",restart_name,OutputID);
        strcat(restart_name,".dat");

        strcpy(roughDensity_name,DirName);
        strcat(roughDensity_name,"/RoughDensity");
        mkdir(roughDensity_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(roughDensity_name,"%s/%d",roughDensity_name,OutputID);
        strcat(roughDensity_name,".dat");

        strcpy(roughDensity1_name,DirName);
        strcat(roughDensity1_name,"/RoughDensity1");
        mkdir(roughDensity1_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(roughDensity1_name,"%s/%d",roughDensity1_name,OutputID);
        strcat(roughDensity1_name,".dat");

        strcpy(roughDensity2_name,DirName);
        strcat(roughDensity2_name,"/RoughDensity2");
        mkdir(roughDensity2_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(roughDensity2_name,"%s/%d",roughDensity2_name,OutputID);
        strcat(roughDensity2_name,".dat");

        strcpy(density_name,DirName);
        strcat(density_name,"/Density");
        mkdir(density_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(density_name,"%s/%d",density_name,OutputID);
        strcat(density_name,".dat");

        strcpy(density1_name,DirName);
        strcat(density1_name,"/Density1");
        mkdir(density1_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(density1_name,"%s/%d",density1_name,OutputID);
        strcat(density1_name,".dat");

        strcpy(density2_name,DirName);
        strcat(density2_name,"/Density2");
        mkdir(density2_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(density2_name,"%s/%d",density2_name,OutputID);
        strcat(density2_name,".dat");

        strcpy(walldensity_name,DirName);
        strcat(walldensity_name,"/WallDensity");
        mkdir(walldensity_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(walldensity_name,"%s/%d",walldensity_name,OutputID);
        strcat(walldensity_name,".dat");

        strcpy(walldensity1_name,DirName);
        strcat(walldensity1_name,"/WallDensity1");
        mkdir(walldensity1_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(walldensity1_name,"%s/%d",walldensity1_name,OutputID);
        strcat(walldensity1_name,".dat");

        strcpy(walldensity2_name,DirName);
        strcat(walldensity2_name,"/WallDensity2");
        mkdir(walldensity2_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(walldensity2_name,"%s/%d",walldensity2_name,OutputID);
        strcat(walldensity2_name,".dat");

        strcpy(roughHeight_name,DirName);
        strcat(roughHeight_name,"/RoughHeight");
        mkdir(roughHeight_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(roughHeight_name,"%s/%d",roughHeight_name,OutputID);
        strcat(roughHeight_name,".dat");

        strcpy(height_name,DirName);
        strcat(height_name,"/Height");
        mkdir(height_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(height_name,"%s/%d",height_name,OutputID);
        strcat(height_name,".dat");

        strcpy(normal_name,DirName);
        strcat(normal_name,"/Normal");
        mkdir(normal_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(normal_name,"%s/%d",normal_name,OutputID);
        strcat(normal_name,".dat");

        strcpy(env_name,DirName);
        strcat(env_name,"/Environment");
        mkdir(env_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(env_name,"%s/%d",env_name,OutputID);
        strcat(env_name,".dat");

        strcpy(aga_name,DirName);
        strcat(aga_name,"/AgarField");
        mkdir(aga_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(aga_name,"%s/%d",aga_name,OutputID);
        strcat(aga_name,".dat");

        strcpy(wal_name,DirName);
        strcat(wal_name,"/WallField");
        mkdir(wal_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(wal_name,"%s/%d",wal_name,OutputID);
        strcat(wal_name,".dat");

	// open files for output
	if (append) Files.cells = fopen(cell_name, "a");	// file for cell statistics output
	else Files.cells = fopen(cell_name, "w");

	if (Files.cells == NULL)
	{
	  std::cout << "Can't open output file." << std::endl;
	  exit(1);
	}

	Files.restart = fopen(restart_name, "w");
	if (Files.restart == NULL) {
	  std::cout << "Can't open restart file." << std::endl;
	  exit(1);
	}

	Files.roughDensity = fopen(roughDensity_name, "w");	// file to store roughDensity of cells
	if (Files.roughDensity == NULL) {
	  std::cout << "Can't open roughDensity file." << std::endl;
	  exit(1);
	}

	Files.roughDensity1 = fopen(roughDensity1_name, "w");	// file to store roughDensity1 of cells
	if (Files.roughDensity1 == NULL) {
	  std::cout << "Can't open roughDensity1 file." << std::endl;
	  exit(1);
	}

	Files.roughDensity2 = fopen(roughDensity2_name, "w");	// file to store roughDensity2 of cells
	if (Files.roughDensity2 == NULL) {
	  std::cout << "Can't open roughDensity2 file." << std::endl;
	  exit(1);
	}

	Files.density = fopen(density_name, "w");	// file to store density of cells
	if (Files.density == NULL) {
	  std::cout << "Can't open density file." << std::endl;
	  exit(1);
	}

	Files.density1 = fopen(density1_name, "w");	// file to store density1 of cells
	if (Files.density1 == NULL) {
	  std::cout << "Can't open density1 file." << std::endl;
	  exit(1);
	}

	Files.density2 = fopen(density2_name, "w");	// file to store density2 of cells
	if (Files.density2 == NULL) {
	  std::cout << "Can't open density2 file." << std::endl;
	  exit(1);
	}

	Files.walldensity = fopen(walldensity_name, "w");	// file to store density of cells
	if (Files.walldensity == NULL) {
		std::cout << "Can't open walldensity file." << std::endl;
	  exit(1);
	}

	Files.walldensity1 = fopen(walldensity1_name, "w");	// file to store density of cells
	if (Files.walldensity1 == NULL) {
		std::cout << "Can't open walldensity1 file." << std::endl;
		 exit(1);
	}

	Files.walldensity2 = fopen(walldensity2_name, "w");	// file to store density of cells
	if (Files.walldensity2 == NULL) {
		std::cout << "Can't open walldensity2 file." << std::endl;
		 exit(1);
	}

        Files.roughheight = fopen(roughHeight_name, "w");	// file to store height of cells
        if (Files.roughheight == NULL) {
            std::cout << "Can't open roughheight file." << std::endl;
            exit(1);
        }
    
	Files.height = fopen(height_name, "w");	// file to store height of cells
	if (Files.height == NULL) {
	  std::cout << "Can't open height file." << std::endl;
	  exit(1);
	}

	Files.normal = fopen(normal_name, "w");	// file to store surface tension forces
	if (Files.normal == NULL) {
	  std::cout << "Can't open surface tension file." << std::endl;
	  exit(1);
	}

	Files.env = fopen(env_name, "w");	// file to store surface tension forces
	if (Files.env == NULL) {
	  std::cout << "Can't open environment file." << std::endl;
	  exit(1);
	}

	Files.aga = fopen(aga_name, "w");	// file to store surface tension forces
	if (Files.aga == NULL) {
	  std::cout << "Can't open agar field file." << std::endl;
	  exit(1);
	}

	Files.wal = fopen(wal_name, "w");	// file to store surface tension forces
	if (Files.wal == NULL) {
	  std::cout << "Can't open wall field file." << std::endl;
	  exit(1);
	}
    

}

void CloseOutputFiles(OutputFiles& Files)
{
        fclose(Files.cells);
        fclose(Files.roughDensity);
        fclose(Files.roughDensity1);
        fclose(Files.roughDensity2);
        fclose(Files.density);
        fclose(Files.density1);
        fclose(Files.density2);
        fclose(Files.walldensity);
        fclose(Files.walldensity1);
        fclose(Files.walldensity2);
        fclose(Files.roughheight);
        fclose(Files.height);
        fclose(Files.env);
        fclose(Files.aga);
        fclose(Files.wal);
        fclose(Files.restart);
        fclose(Files.normal);
}

int AddDropCells(Cell* cells, double L_divide, double radius, UniformGrid& Grid, Inputs& Ini)
{
        double PI = 4.0*atan(1.0);
        double dz = 0.0;
        int icell = 0;

        double L, thetaPos;
        double thetaDir, radiusPos;
        DoubleCoord v, va, p, q, cm, c1, c2;
        bool CheckOverlap = true;
        int RegenCellMax = 10000;
        double dist;

        // Construct my random seed, based on time
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        std::uniform_real_distribution<double> distribution (0.0,1.0);

        while (icell<Ini.ColonySize)
        {
                int RegenCellCount = 0;
                CheckOverlap = true;
                while (CheckOverlap==true)
                {
                        L                   = L_divide*(max(0.3,0.95*distribution(generator)));
                        cells[icell].Length = L;
                        cells[icell].Radius = radius;
                        radiusPos           = distribution(generator)*Ini.ColonyRadius;
                        thetaPos            = 2.0*PI*drand48();
                        thetaDir            = PI*drand48();

                        cm = DoubleCoord(radiusPos*cos(thetaPos), radiusPos*sin(thetaPos), radius+dz);
                        p  = DoubleCoord( 0.5*L*cos(thetaDir) + cm.x ,  0.5*L*sin(thetaDir) + cm.y, cm.z);
                        q  = DoubleCoord(-0.5*L*cos(thetaDir) + cm.x , -0.5*L*sin(thetaDir) + cm.y, cm.z);
                        v  = DoubleCoord(0,0,0);
                        va = DoubleCoord(0,0,0);
                        cells[icell].Radius          = radius;
                        cells[icell].Position.p      = p;
                        cells[icell].Position.q      = q;
                        cells[icell].Position.time_p = 0;
                        cells[icell].Position.time_q = 0;
                        cells[icell].Position.age_p  = 0;
                        cells[icell].Position.age_q  = 0;
                        cells[icell].Velocity        = v;
                        cells[icell].AngularVelocity = va;
                        cells[icell].Ancestor        = icell+1;

                        cells[icell].GrowthRate      = maxGrowthRate;
                        cells[icell].AgeCell         = 0.0;
                        cells[icell].Type            = 1;
                        cells[icell].Ldiv            = L_divide;

                        Grid.Add(icell, Grid.GetAddress(cm));

                        int icheck=0;
                        while (icheck<icell)
                        {
                        min_distance(cells[icell],cells[icheck],dist,c1,c2);
                                if (dist<(cells[icell].Radius+cells[icheck].Radius))
                                {
                                        CheckOverlap=true;
                                        //std::cout << "Cells overlap!" << std::endl;
                                        RegenCellCount++;
                                        break;
                                }
                                icheck = icheck + 1;
                        }
                        if (icheck==icell) {CheckOverlap=false;}
                        if (RegenCellCount==RegenCellMax)
                        {
                                std::cout << "Unable to generate initial cells!" << std::endl;
                                exit(0);
                        }
                }
                icell = icell + 1;
        }

        std::cout << "Initial cells drop-like placed" << std::endl;
        return icell;
}

int AddFirstCells(Cell* cells, double L_divide, double radius, UniformGrid& Grid, Inputs& Ini, const int caseTest)
{
        int icell = 0;
        icell     = AddDropCells(cells, L_divide, radius, Grid, Ini);
        t0        = 0;
        return icell;
}

int LoadCells(char* fname, Cell* cells, UniformGrid& Grid, double& t, double& dt)
{

	std::cout << "Reading cells from " << fname << std::endl;

	FILE* FID = fopen(fname, "r");
	if (FID == NULL) {
	  std::cout << "Can't open restart file." << std::endl;
	  exit(1);
	}

	// obtain file size:
	fseek (FID, 0, SEEK_END);
	int fsize = ftell (FID);
	rewind (FID);

	// read time
	fread(&t, sizeof(double), 1, FID);
	fread(&dt, sizeof(double), 1, FID);

	std::cout << "(source) t = " << std::setprecision(4) << t << ", (source) dt = " \
	          << std::setprecision(4) << dt << std::endl;

	t0 = t;
	// read cells
	int cell_count = (fsize-2*sizeof(double))/sizeof(Cell);

	fread (cells, sizeof(Cell), cell_count, FID);
	std::cout << "Read " << cell_count << " cells" << std::endl;

        DoubleCoord v, va;
        v = DoubleCoord(0,0,0);
        va = DoubleCoord(0,0,0);
	
	for (int icell = 0; icell<cell_count; icell++)
	{
		Grid.Add(icell, Grid.GetAddress(average(cells[icell].Position)));
		cells[icell].Velocity = v;
		cells[icell].AngularVelocity = va;
		cells[icell].GrowthRate = maxGrowthRate;
		cells[icell].Ldiv = L_divide;
	}
	std::cout << "Added to grid " << std::endl;

	fclose(FID);

	return cell_count;
}

void SaveCells(FILE* FID, Cell* cells, int N_cells, double t, double dt)
{
	// save cell information
	rewind(FID);

	int size_written = 0;

	size_written = fwrite(&t, sizeof(double), 1, FID);
	//MyAssert(size_written>0,"Could not write restart file");

	fwrite(&dt, sizeof(double), 1, FID);
	fwrite(cells, sizeof(Cell), N_cells, FID );
	fflush(FID);
}

Inputs ReadParameters(char* fname)
{
	FILE* FID = fopen(fname, "r");
	if (FID == NULL) 
	{
	        std::cout << "Can't open parameter file." << std::endl;
	        exit(1);
	}
	
	char* data_string;
	char var_name[100];
	char var_value[100];

	int fileLen = GetFileLen(FID);
	char* buffer = (char*) malloc(fileLen+1);
	fread(buffer, fileLen, 1, FID);
	buffer[fileLen] = 0;

	Inputs IniConditions;

	std::cout << "Reading: " << fname << std::endl;
	
	while(data_string = GetNextString(buffer))
	{

	//while (fscanf(FID, "%s %f \r", var_name, var_value) != NULL)
	//while (fgets (data_string , 100 , FID) != NULL)
	//{
		sscanf(data_string, "%s %s", var_name, var_value);
		std::cout << "\t" << var_name << "\t" << var_value << std::endl;

		if (strcmp(var_name,"Radius")==0)
			cellRadius = atof(var_value);
		else if (strcmp(var_name,"L_divide")==0)
			L_divide = atof(var_value);
		else if (strcmp(var_name,"k_cc")==0)
			k_cc = atof(var_value);
		else if (strcmp(var_name,"k_wc")==0)
			k_wc = atof(var_value);
		else if (strcmp(var_name,"var_L")==0)
			varL = atof(var_value);
		else if (strcmp(var_name,"var_angle")==0)
			varAngle = atof(var_value);
                else if (strcmp(var_name,"var_pos")==0)
                        var_pos = atof(var_value);
		else if (strcmp(var_name,"Viscosity")==0)
			viscosity = atof(var_value);
		else if (strcmp(var_name,"Growth_Rate")==0)
			maxGrowthRate = atof(var_value);
		else if (strcmp(var_name,"Wall_Rough")==0)
			wall_rough = atof(var_value);
		else if (strcmp(var_name,"Gamma_t")==0)
			gamma_t = atof(var_value);
		else if (strcmp(var_name,"Gamma_n")==0)
			gamma_n = atof(var_value);
		else if (strcmp(var_name,"Wall_Mu")==0)
			wall_mu = atof(var_value);
		else if (strcmp(var_name,"Cell_Mu")==0)
			cell_mu = atof(var_value);
		else if (strcmp(var_name,"Density_Threshold")==0)
			density_threshold = atof(var_value);
		else if (strcmp(var_name,"Surface_Tension")==0)
			tension = atof(var_value);
		else if (strcmp(var_name,"t_max")==0)
			t_max = atof(var_value); 
		else if (strcmp(var_name,"dt")==0)
			initial_dt = atof(var_value);	
		else if (strcmp(var_name,"Box_x")==0)
			BoxX = atoi(var_value);
		else if (strcmp(var_name,"Box_y")==0)
			BoxY = atoi(var_value);
		else if (strcmp(var_name,"Box_z")==0)
			BoxZ = atoi(var_value);
		else if (strcmp(var_name,"Box_z_agar")==0)
			BoxZAgar = atoi(var_value);
		else if (strcmp(var_name,"Box_Dim")==0)
		{
			BoxX = atoi(var_value);
			BoxY = BoxX;
		}
		else if (strcmp(var_name,"BoxLength")==0)
			BoxLength = atof(var_value); 
                else if (strcmp(var_name,"maxLevels")==0)
                        maxLevels = atoi(var_value);
                else if (strcmp(var_name,"refinementGridHeight")==0)
                        refinementGridHeight = atoi(var_value);
		else if (strcmp(var_name,"Output_Time")==0)
			OutputTime = atof(var_value);
		else if (strcmp(var_name,"Update_Time")==0)
			UpdateTime = atof(var_value);
		else if (strcmp(var_name,"Tortuosity")==0)
			Tortuosity = atof(var_value);
		else if (strcmp(var_name,"KC")==0)
			KC = atof(var_value);
		else if (strcmp(var_name,"C_rate")==0)
			C_rate = atof(var_value);
		else if (strcmp(var_name,"Diff_Colony")==0)
			DiffColony = atof(var_value);
		else if (strcmp(var_name,"Diff_Agar")==0)
			DiffAgar = atof(var_value);
		else if (strcmp(var_name,"maxCarbon")==0)
			maxCarbon = atof(var_value);
		else if (strcmp(var_name,"Cdt")==0)
                        Cdt = atof(var_value);
		else if (strcmp(var_name,"ConvCrit")==0)
			ConvCrit = atof(var_value);
		else if (strcmp(var_name,"minIter")==0)
			minIter = atof(var_value);
		else if (strcmp(var_name,"maxIter")==0)
			maxIter = atof(var_value);
		else if (strcmp(var_name,"InterfaceCondition")==0)
			InterfaceCondition = atof(var_value);
		else if (strcmp(var_name,"NutrientGSI")==0)
			NutrientGSI = (bool)atoi(var_value);
		else if (strcmp(var_name,"Rc")==0)
			Rc = atof(var_value);
		else if (strcmp(var_name,"Delta_H")==0)
			DH = atof(var_value);
		else if (strcmp(var_name,"MaintenanceRate")==0)
			Maintenance_rate = atof(var_value);
		else if (strcmp(var_name,"FilterLen")==0)
			FilterLen = atoi(var_value);
		
		// Initial condition stuff
		else if (strcmp(var_name,"IniColonyRadius")==0)
			IniConditions.ColonyRadius = atof(var_value);
		else if (strcmp(var_name,"IniColonySize")==0)
			IniConditions.ColonySize = atoi(var_value);
		else if (strcmp(var_name,"NumColonies")==0)
			IniConditions.ColonyNumber = atoi(var_value);
		else if (strcmp(var_name,"ColonySeparation")==0)
			IniConditions.ColonySeparation = atof(var_value);
                else if (strcmp(var_name,"MaxCells")==0)
                        maxCells = atoi(var_value);

                // Phage stuff added here
                else if (strcmp(var_name,"latent_period")==0)
                        latent_period = atof(var_value);
                else if (strcmp(var_name,"colony_size")==0)
                        colony_size = atoi(var_value);
                else if (strcmp(var_name,"colony_time")==0)
                        colony_time = atof(var_value);
                else if (strcmp(var_name,"ProlePhage")==0)
                        ProlePhage = atoi(var_value);
                else if (strcmp(var_name,"Abortive_frac")==0)
                        Abortive_frac = atof(var_value);
                else if (strcmp(var_name,"Infect_frac")==0)
                        Infect_frac = atof(var_value);
                else if (strcmp(var_name,"Adsorption_rate")==0)
                        Adsorption_rate = atof(var_value);

		else
		{
			std::cout << "Unknown parameter: " << var_name << std::endl;
			fflush(stdout);
			assert(false);
			exit(-1);
		}
	}
	std::cout << "Done Reading: " << fname << std::endl << std::endl;
	fclose(FID);
	return IniConditions;
}

int GetFileLen(FILE* myFile)
{
	fseek (myFile, 0, SEEK_END);
	int size = ftell(myFile);
	fseek(myFile, 0, SEEK_SET);
	return size;
}

char* GetNextString(char*& buffer)
{
        char* out = buffer;
        if (!*buffer) return NULL; // return on empty string
        while(! (*buffer == 0x0A || *buffer == 0x0D || *buffer == 0x00) ) // 0x0A and 0x0D
                buffer++; // skip forward until we find the start of the next line (10/13/0)
        if (*buffer) *buffer++ = 0; // if we ended on 10/13 end the string and move to the next char
        if(*buffer == 0x0A) buffer++;  // on windows skip the 10 after the 13

        return out;
}

void Output(FILE* FID, int ID, double t, const Cell& cell, const Tensor T)
{

	int Type;
	if(cell.Short)
	{
	        Type = 2;
	}
	else
	{
	        Type = 1;
	}
	
	if(cell.Shrink and cell.Resistant == false)
	{
	        Type = 3;
	}
	
	if(cell.Resistant)
	{
	        Type = 4;
	}
	if(cell.PhageCell > 0)
	{
	        Type = 66;
	}
	if(cell.PhageCell > 66 and cell.Resistant == false)
	{
	        Type = 77;
	}
	
	{
		fprintf(FID,"%d %d %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E\n",
		ID, Type, cell.PhageCell, \
		cell.Position.p.x, cell.Position.p.y, cell.Position.p.z, \
		cell.Position.q.x, cell.Position.q.y, cell.Position.q.z, \
		cell.Length, T.xx, T.yy, T.zz, \
		cell.Velocity.x, cell.Velocity.y, cell.Velocity.z, \
		cell.GrowthRate,cell.AgeCell,cell.Pressure,cell.Strain);
	}
}

void Output(FILE* FID, int ID, double t, const Cell& cell, const DoubleCoord F)
{

	int Type;
	if(cell.Short)
	{
	        Type = 2;
	}
	else
	{
	        Type = 1;
	}
	
	if(cell.Shrink and cell.Resistant == false)
	{
	        Type = 3;
	}
	
	if(cell.Resistant)
	{
	        Type = 4;
	}
	if(cell.PhageCell > 0)
	{
	        Type = 66;
	}
	if(cell.PhageCell > 66 and cell.Resistant == false)
	{
	        Type = 77;
	}
	
	{
	        fprintf(FID,"%d %d %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E %4.4E\n",
                ID, Type, cell.PhageCell, \
                cell.Position.p.x, cell.Position.p.y, cell.Position.p.z, \
                cell.Position.q.x, cell.Position.q.y, cell.Position.q.z, \
                cell.Length, F.x, F.y, F.z, \
                cell.Velocity.x, cell.Velocity.y, cell.Velocity.z, \
                cell.GrowthRate,cell.AgeCell,cell.Pressure,cell.Strain);
	}
}



