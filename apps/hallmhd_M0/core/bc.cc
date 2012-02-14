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
    int a;

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
        // axis boundary condition for a radial z-pinch
          
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; cc<=ncoeffs; cc++)
            {
                if (cc%2!=0) a=1;
                else a=-1;

                // scalars are copied over to ghost cells
                 q(1-ibc,1,cc) = a*q(ibc,1,cc);
                 q(1-ibc,5,cc) = a*q(ibc,5,cc);
                 q(1-ibc,6,cc) = a*q(ibc,6,cc);
                 q(1-ibc,10,cc) = a*q(ibc,10,cc);

                // radial terms are flipped in sign
                q(1-ibc,2,cc) = -a*q(ibc,2,cc);
                q(1-ibc,7,cc) = -a*q(ibc,7,cc);
                q(1-ibc,11,cc) = -a*q(ibc,11,cc);
                q(1-ibc,14,cc) = -a*q(ibc,14,cc);

                // phi terms are flipped in sign
                q(1-ibc,3,cc) = -a*q(ibc,3,cc);
                q(1-ibc,8,cc) = -a*q(ibc,8,cc);
                q(1-ibc,12,cc) = -a*q(ibc,12,cc);
                q(1-ibc,15,cc) = -a*q(ibc,15,cc);

                // z terms are copied across axis
                 q(1-ibc,4,cc) = a*q(ibc,4,cc);
                 q(1-ibc,9,cc) = a*q(ibc,9,cc);
                 q(1-ibc,13,cc) = a*q(ibc,13,cc);
                 q(1-ibc,16,cc) = a*q(ibc,16,cc);

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
        // conducting wall boundary condition for a radial z-pinch
        // reflect all fluid and EM variables
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                {
                    if (cc%2!=0) a=1;
                    else a=-1;
                    q(mx+ibc,m,cc) = a*q(mx+1-ibc,m,cc);
                }

        // negate fluid normal velocities
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; cc<=ncoeffs; cc++)
            {
                if (cc%2!=0) a=1;
                else a=-1;
                q(mx+ibc,2,cc) = -a*q(mx+1-ibc,2,cc);
                q(mx+ibc,7,cc) = -a*q(mx+1-ibc,7,cc);
            }

        // set correct signs for EM fields
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; cc<=ncoeffs; cc++)
            {
                if (cc%2!=0) a=1;
                else a=-1;

                // negate Ey,Ez,Bx
                q(mx+ibc,12,cc) = -a*q(mx+1-ibc,12,cc);
                q(mx+ibc,13,cc) = -a*q(mx+1-ibc,13,cc);
                q(mx+ibc,14,cc) = -a*q(mx+1-ibc,14,cc);
            }
    }

}

