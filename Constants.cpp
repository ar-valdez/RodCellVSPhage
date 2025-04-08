#include <math.h>

// Mechanical parameters (Table S4 Warren2019.pdf)
double cell_mu                = 1.0e-01;
double wall_mu                = 8.0e-01;
double k_cc                   = 3.0e+04;// elastic constant for cell-cell interactions (atm?)
double k_wc                   = 3.0e+04;// elastic constant between cells and wall
double gamma_t                = 1.0e+04;// tangential dissipation rate
double gamma_n                = 5.0e+02;// normal dissipation rate 100
double tension                = 1.5e+02;// surface tension
double viscosity              = 1.00160e-03;// water viscosity https://en.wikipedia.org/wiki/List_of_viscosities


// E. Coli parameters (Table S3 Warren2019.pdf)
double cellRadius             = 5.0e-01;// radius of each cell cap (constant for now, microns)
double L_divide               = 3.0;// 3 length when cells divide (microns)
double varL                   = 0.025;// max variation in the length of daughter cells
//double varL                   = 9.419e-01;// Fig 2b, Gangan2017.pdf
double varAngle               = 5.0e-04;// max variation in the orientation of daughter cells
double maxGrowthRate          = 8.8e-01;// Growth rate Arabinose Table S1
double density_threshold      = 0.68;//1.02e-01;// Mya set this value to 0.6 previously.

double var_pos                = 0.0;
double wall_rough             = 0.0;
double DH                     = -1.0e-02;// determines how tightly the surface tension holds the cells (+ value is low agar concentration, - value is high agar concentration)

// Nutrient parameters (Table S2 Warren2019.pdf)
double Tortuosity             = sqrt(2.0);
double KC                     = 2.0e+01;// Monods' Nutrient capacity constant (K_S Table S2)
double C_rate                 = 3.7e-5*Tortuosity*5.0;
double maxCarbon              = 1.33e+01;// Initial concentration of Carbon Arabinose case
double DiffColony             = 1.5e-01;
double DiffAgar               = 1.0e+00;

double Maintenance_rate       = 0.0;
double Rc                     = 3.0;

// Temporal parameters 
double t0                     = 0.0;// initial time
double t_max                  = 20.0;// total simulation time
double initial_dt             = 2.0e-05;// initial time step
double OutputTime             = 1.0e-01;// how often to output
double UpdateTime             = 3.0e-01;// how often to solve nutrients

// Geometry parameters
int BoxX                      = 300;//128;
int BoxY                      = 300;//128;
int BoxZ                      = 10;//20;
int BoxZAgar                  = 10;//16;
int maxLevels                 = 4;
double BoxLength              = 2.0*(L_divide+2*cellRadius);
int FilterLen                 = 5;

// Other parameters
int maxCells                  = 100000;
double Cdt                    = 0.1;
double ConvCrit               = 5.0e-02;
int minIter                   = 10;
int maxIter                   = 10000;
int InterfaceCondition        = 1;//1: continuous $\partial C/partial n$
//int InterfaceCondition        = 2;//2: flux continuity with qC
//int InterfaceCondition        = 3;//3: continuous $\partial C/\partial t$
bool NutrientGSI              = 0;
int refinementGridHeight      = 2;

// Phage parameters
double latent_period          = 0.5;  // time to experience burst (hours)
int colony_size               = 1000; // Min size of colony to start infection
double colony_time            = 100;   // colony time before phage infection
int ProlePhage                = 300;  // phages generated after burst event
double Abortive_frac          = 0.0;  // Fraction of abortive cells
double Infect_frac            = 3.0;  // Fraction of infected cells
double Adsorption_rate        = 3.0e-09; // Adsorption Rate 
// directory name
char DirName[500]             = "";

// Phage diffusion parameters
double PI                     = 4.0*atan(1);
double eta                    = 5.0e+04; // Pa. s, Viscosity at  298.15 K
double Kb                     = 1.380649e-26; // Boltzman cts. @ SI
double temp                   = 298.15;; //Temperature Kelvin
double ph_length              = 200e-09; // Length of a Phage
double Dph                    = Kb * temp / (6.0*PI*eta*ph_length); //Einstein diffusion

// Case handling
int caseTest                  = 1; // No phage
//int caseTest                  = 2; // PacMan
//int caseTest                  = 3; // Border Cells
//int caseTest                  = 4; // Plaque Assay
//int caseTest                  = 5; // Double PacMan
//int caseTest                  = 6; // Bug case
//int caseTest                  = 7; // SnowFlake
//int caseTest                  = 8; // Mushroom
//int caseTest                  = 9; // HalfMoon
//int caseTest                  = 10; // 4 cells at axis
//int caseTest                  = 11; // Cells random top
//int caseTest                  = 12; // Phage random uniform
//int caseTest                  = 13; // Center Cells

//int caseTest                  = 20; // Kaz killing case
//int caseTest                  = 25; // Kaz longFriction case
//int caseTest                  = 30; // Exponential
//int caseTest                  = 40; // Hertzian
//int caseTest                  = 35; // MixForce

// Chamber values
double Vwall                  = 100.0;
double Hwall                  = 0.00;
double Hgtflush               = 40.0;


