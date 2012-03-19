#include <iostream>
#include <fstream>
#include <math.h>

#include "farray.h"
#include "miniwarpx.h"


static const int LEFT = 0;
static const int RIGHT = 1;

/**
   Standard boundary routine for WARPX. Choices which work for all
   equation systems and algorithms are bc_periodic and bc_copy. Other
   choices like bc_wall and bc_axis are highly equation and algorithm
   specific and should be handled by users by making a copy of this
   file in their own directories and then modifing the code in this
   function.

*/
void 
bc(const Run_Data& rd, FArray<double>& q)
{
    const int *bc_type = rd.bc_type;
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mbc = rd.mbc;
    int ncoeffs = rd.ncoeffs;

    //
    // Apply boundary condition to the left boundary
    //
    if(bc_type[LEFT] == bc_periodic)
        // periodic boundary conditions
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(1-ibc,m,cc) = q(mx-ibc+1,m,cc);

    else if(bc_type[LEFT] == bc_copy)
        // zero order extrapolation
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(1-ibc,m,cc) = q(1,m,cc);
    
    else if(bc_type[LEFT] == bc_custom)
    {
        double t = rd.tcurrent;
        double t0 = 1.0;
        double D = 0.5;

        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(1-ibc,m,cc) = 0.0;

        for(int ibc=1; ibc<=mbc; ibc++)
        {
            // simulate an incoming wave
            
            q(1-ibc,2,1) = exp(-pow(t-t0,2)/(D*D))*sin(20*t); // Ey
            q(1-ibc,3,1) = exp(-pow(t-t0,2)/(D*D))*cos(20*t); // Ez

            q(1-ibc,5,1) = -q(1-ibc,3,1); // By = Ez
            q(1-ibc,6,1) = q(1-ibc,2,1); // Bz = Ey
        }
    }

    //
    // Apply boundary condition to the right boundary
    //
    if(bc_type[RIGHT] == bc_periodic)
        // periodic boundary conditions
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(mx+ibc,m,cc) = q(ibc,m,cc);
    
    else if(bc_type[RIGHT] == bc_copy)
        // zero order extrapolation
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(mx+ibc,m,cc) = q(mx,m,cc);

    else if(bc_type[RIGHT] == bc_custom)
    {
        // do your own thing
    }
}

