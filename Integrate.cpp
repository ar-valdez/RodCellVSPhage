#include <math.h>
#include <iostream>
//#include <cmath.h>

#include "Cell.h"
#include "Forces.h"
#include "Integrate.h"
#include "tools.h"

// Perform the numerical integration
// The motion Eqs. are:
// \eta\,\dt{x}      = F,
// \eta\,\dt{\theta} = r x F = T,

void integrate(const double dt, int cellID, const Cell* old_cells, Cell* new_cells, \
               const int* neighbours, DoubleArray2D& Height, CoordArray2D& Normal, \
               UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Wall, \
               const int caseTest)
{
        DoubleCoord F, T;
        double press  = 0.0;
        
        // Heuns method
        
        // cell we're computing forces on
        Cell temp_tildecell;
        DoubleCoord tilde_F, tilde_T;
        
        // Get forces at source
        sum_forces(old_cells[cellID], old_cells, neighbours, F, T, Height, Normal, Grid, XYAddress, Wall, caseTest, press);

        // Solve position for tilde step
        getPositions(dt, F, T, old_cells[cellID], temp_tildecell);

        // Get forces at tilde position
        sum_forces(temp_tildecell, old_cells, neighbours, tilde_F, tilde_T, Height, Normal, Grid, XYAddress, Wall, caseTest, press);
        
        // Solve new positions using tilde Forces and original forces
        F = sum(F,tilde_F);
        T = sum(T,tilde_T);
        
        // Initial Guess (From Heuns)
        getPositions(0.5*dt, F, T, old_cells[cellID],new_cells[cellID]);
        new_cells[cellID].Pressure = press;
        
}

// Here I Update positions using new overdamped equation for translation and rotation
// Velocities are no longer needed.
void getPositions(const double dt, const DoubleCoord& F, const DoubleCoord& T, const Cell& old_cell, Cell& new_cell)
{
        // Get a target cell to work in
        new_cell = old_cell;

        double eta = 1.0e+06 * viscosity;
        double pc  = (old_cell.Length + old_cell.Radius) / old_cell.Radius;
        double ct  = 0.312 + (0.565/pc) - (0.1/pow(pc,2));
        double cr  = -0.662 + (0.917/pc) - (0.05/pow(pc,2));
        double PI  = 4.0*atan(1.0);
        
        double eta_t = (3.0*PI*eta*(old_cell.Length + old_cell.Radius)) / (log(pc) + ct);
        double eta_r = (PI*eta/pow(old_cell.Length,3)) / (3*log(pc)+3*cr);
        
        DoubleCoord Feval = F;
        DoubleCoord Teval = T;
        
        // Get new translational speed
        new_cell.Velocity        = scale(Feval,(1.0/eta_t)); // Lev's advice + Igor's Newton advice

        // Get new rotational speed
        new_cell.AngularVelocity = scale(Teval,(1.0/eta_t)); // Lev's advice + Igor's Newton advice + eta_t test
        
        // Work with centre of mass (Rigid motion)
        DoubleCoord cm = average(old_cell.Position);
        
        // Get relative locations of (p,q)
        DoubleCoord rp = diff(old_cell.Position.p,cm);
        DoubleCoord rq = diff(old_cell.Position.q,cm);

        // Get Angular-Translation speed on (p,q)
        DoubleCoord vrp = cross(new_cell.AngularVelocity,rp);
        DoubleCoord vrq = cross(new_cell.AngularVelocity,rq);

        // Get total velocity
        DoubleCoord vp = sum(new_cell.Velocity,vrp);
        DoubleCoord vq = sum(new_cell.Velocity,vrq);
        
        // Integrate to find new positions
        new_cell.Position.p = sum(old_cell.Position.p, scale(vp,dt));
        new_cell.Position.q = sum(old_cell.Position.q, scale(vq,dt));

        //new_cell.Position.p.z = new_cell.Radius;
        //new_cell.Position.q.z = new_cell.Radius;
        
        // Fixings to avoid Cell extension to length

        // New cm
        DoubleCoord cm_new = average(new_cell.Position);
        
        // Get relative locations of (p,q)
        rp = diff(new_cell.Position.p,cm_new);
        rq = diff(new_cell.Position.q,cm_new);
        
        //new cell length
        DoubleCoord uv = diff(new_cell.Position.p,new_cell.Position.q);
        
        double factor = old_cell.Length / sqrt(dot(uv,uv)); 
        new_cell.Position.p = sum(cm_new, scale(rp,factor));
        new_cell.Position.q = sum(cm_new, scale(rq,factor));
       
}




