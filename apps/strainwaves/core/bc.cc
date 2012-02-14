#include <iostream>
#include <fstream>
#include <math.h>
#include "strainwave.h"
#include "farray.h"
#include "miniwarpx.h"
using namespace std;

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
    Strainwave_Vars *sv = (Strainwave_Vars*) rd.mvar;
    double u;                      // sets boundary condition of velocity
    double t = rd.tcurrent;        // current time of simulation
    double pi = 3.14159;

    if(bc_type[LEFT] == bc_periodic)
        // periodic boundary conditions
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(1-ibc,m,cc) = q(mx+ibc-1,m,cc);

    else if(bc_type[LEFT] == bc_copy)
        // zero order extrapolation
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(1-ibc,m,cc) = q(1,m,cc);
    
    else if(bc_type[LEFT] == bc_custom)
    {
        if (t <= 20)
        {
            u = -sv->ubar*(1+cos(pi*(t-10)/10));
        }
        else
        {
            u=0;
        }
        // inflow boundary condition on the left
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; cc<=ncoeffs; cc++)
            {
                q(2-ibc,1,cc) = 0.0;
                q(2-ibc,2,cc) = sv->rho(2-ibc)*u;
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
        // outflow boundary condition on the right
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; cc<=ncoeffs; cc++)
            {
            }
    }
}

