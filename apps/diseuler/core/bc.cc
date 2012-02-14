#include <iostream>
#include <fstream>
#include <math.h>
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
        // axis boundary condition
          
        double a;
        for(int ibc=1; ibc<=mbc; ibc++)
        {
            for(int cc=1; cc<=ncoeffs; cc++)
            {
                a = pow(-1.,cc-1);

                // scalars are copied over to ghost cells
                 q(1-ibc,1,cc) = a*q(ibc,1,cc);
                 q(1-ibc,5,cc) = a*q(ibc,5,cc);

                // radial terms are flipped in sign
                q(1-ibc,2,cc) = -a*q(ibc,2,cc);

                // phi terms are flipped in sign
                q(1-ibc,3,cc) = -a*q(ibc,3,cc);

                // z terms are copied across axis
                 q(1-ibc,4,cc) = a*q(ibc,4,cc);
            }
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

    }

}

