#include <math.h>
#include <stdlib.h>
#include <float.h>

#include <chrono>
#include <random>

#include "Array.h"
#include "Cell.h"
#include "Constants.h"
#include "Forces.h"
#include "grow.h"
#include "tools.h"
#include "UniformGrid.h"
#include "Nutrients.h"
#include "Phages.h"

using namespace std;

// grow the cell
void grow(double dt, Cell& cell, EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid, bool Nutrients)
{
        if(cell.PhageCell == 0)
        {
                // Compute the real Growth rate (convert liquid to solid)
                double sigma_Ecoli   = log( (2.0*cell.Ldiv) / (cell.Ldiv - 2.0*cell.Radius) ) / log(2);
                double PenaltyGrowth = 1.0;
                if(cell.Abortive or cell.Resistant)
                {
                        PenaltyGrowth = 0.9; //The cost of being better
                }
                
                if(cell.Shrink)
                {
                        PenaltyGrowth = 0.0;
                }
                
                double Growth_rate = 1.0;
                
                if(Nutrients)
                {
                        // Affect Growth rate to nutrient availability
                        // Get position in uniform grid to access correct index for height
                        DoubleCoord cm = average(cell.Position);	// center of mass
                        IntCoord XYAddress = Grid.GetXY(Grid.GetAddress(cm));

                        // Look up the growth rate in the environment array
                        Growth_rate = Env.Get(XYAddress).GrowthRate;

                        if (XYAddress.z==0)
                        {
                                double lambda=cm.z/BoxLength-floor(cm.z/BoxLength);
                                double y0, y1;
                                double Cgr;
                                y0=(Wal[0]->Get(XYAddress.x,XYAddress.y).CarbonAgar+Wal[0]->Get(XYAddress.x-1,XYAddress.y).CarbonAgar+Wal[0]->Get(XYAddress.x,XYAddress.y-1).CarbonAgar+Wal[0]->Get(XYAddress.x-1,XYAddress.y-1).CarbonAgar)/4;
                                y1=(Env.Get(XYAddress.x,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x,XYAddress.y-1,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y-1,XYAddress.z).Carbon)/4;
                                Cgr = y1*lambda+y0*(1-lambda);
                                Cgr = Cgr/(Cgr+KC);
                                Growth_rate = max(0.0,Cgr-Maintenance_rate/C_rate);
                        }
                        else if (XYAddress.z>0)
                        {
                                double lambda=cm.z/BoxLength-floor(cm.z/BoxLength);
                                double y0, y1;
                                double Cgr;
                                y0=(Env.Get(XYAddress.x,XYAddress.y,XYAddress.z-1).Carbon+Env.Get(XYAddress.x-1,XYAddress.y,XYAddress.z-1).Carbon+Env.Get(XYAddress.x,XYAddress.y-1,XYAddress.z-1).Carbon+Env.Get(XYAddress.x-1,XYAddress.y-1,XYAddress.z-1).Carbon)/4;
                                y1=(Env.Get(XYAddress.x,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y,XYAddress.z).Carbon+Env.Get(XYAddress.x,XYAddress.y-1,XYAddress.z).Carbon+Env.Get(XYAddress.x-1,XYAddress.y-1,XYAddress.z).Carbon)/4;
                                Cgr = y1*lambda+y0*(1-lambda);
                                Cgr = Cgr/(Cgr+KC);
                                Growth_rate = max(0.0,Cgr-Maintenance_rate/C_rate);
                        }

                }
                
                // Finalize and bring the Growth rate back
                Growth_rate = Growth_rate * maxGrowthRate * sigma_Ecoli * PenaltyGrowth;
                
                // Net growth length
                double dL       = cell.Length*dt*Growth_rate;
                
                DoubleCoord v   = diff(cell.Position.q, cell.Position.p);         // vector along segment

                // growth direction
                DoubleCoord dv  = scale(v,dL*0.5/cell.Length);

                cell.Length     = cell.Length + dL; // increase the length of the cell

                cell.Strain     = 100.0*dL / (cell.Ldiv + dL); // Strain at time stamp
                
                cell.Position.p = diff(cell.Position.p,dv);
                cell.Position.q = sum(cell.Position.q,dv);
                cell.GrowthRate = Growth_rate;
                if (cell.Ancestor==0)
                        std::cout << "wrong growth! Ancestor is zero!" << std::endl;
        }
}        


// takes in a mother cell and returns mother and daughter after division
void divide(Cell& mother, Cell& daughter, double t)
{
        if(mother.PhageCell == 0)
        {
                
                // Construct my random seed, based on time
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::default_random_engine generator (seed);
                std::uniform_real_distribution<double> distribution (0.05,0.2);
                double dl      = distribution(generator);               // random part of length after division
                
                DoubleCoord uv = diff(mother.Position.p,mother.Position.q);// vector in between two coordinates of the cell
                double L       = mother.Length;// length of original cell
                uv             = scale(uv,1/L);// scale unit vector by length to get unit vector

                DoubleCoord midpoint = average(mother.Position);
                
                // Now start the division process
                if(mother.Ancestor==0)
                        std::cout << "wrong mother! Ancestor is zero!" << std::endl;
                daughter.Ancestor = mother.Ancestor;
                
                // Division is not deterministic nor equivalent for each new cell
                
                if(mother.Short)
                {
                        mother.Length           = 0.5*L-mother.Radius; // Default criteria Warren2019.pdf
                        daughter.Length         = 0.5*L-mother.Radius;
                }
                else
                {
                        mother.Length           = 0.5*L-mother.Radius+dl; // Default criteria Warren2019.pdf
                        daughter.Length         = 0.5*L-mother.Radius-dl;
                }

                // Update atributes
                daughter.Type           = mother.Type;
                daughter.Abortive       = mother.Abortive;
                daughter.Shrink         = mother.Shrink;
                daughter.Resistant      = mother.Resistant;
                daughter.Short          = mother.Short;
                daughter.PhageCell      = mother.PhageCell;
                daughter.Ldiv           = mother.Ldiv;
                daughter.Radius         = mother.Radius;
                daughter.GrowthRate     = mother.GrowthRate;
                daughter.Velocity       = mother.Velocity;
                daughter.Nb             = mother.Nb;
                daughter.Pressure       = mother.Pressure;
                
                mother.AgeCell          = 0.0;
                daughter.AgeCell        = 0.0;
                
                // update endpoint P,Q for both cells
                
                // Invariant Points
                // p in first  cell is the same as mother cell
                mother.Position.p   = mother.Position.p;
                
                // q in second daughter cell is the same as in the mother cell
                daughter.Position.q = mother.Position.q;
                
                // New points
                // Finish the first cell
                mother.Position.q   = diff(mother.Position.p,scale(uv,mother.Length));
                
                // Finish the second cell
                daughter.Position.p = diff(mother.Position.q, scale(uv,mother.Radius+daughter.Radius));
                //daughter.Position.q = diff(daughter.Position.p, scale(uv,daughter.Length));

                // set new angular velocities

                daughter.AngularVelocity = DoubleCoord(0,0,0);
                mother.AngularVelocity   = DoubleCoord(0,0,0);

                // q in first daughter cell
                daughter.Position.time_q = mother.Position.time_q;
                mother.Position.time_q   = t;
                daughter.Position.time_p = t;
                mother.Position.age_p++;
                daughter.Position.age_q  = mother.Position.age_q+1;
                //mother.Position.time_q = 0;
                //daughter.Position.time_p = 0;
        }
}

