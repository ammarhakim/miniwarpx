#include <iostream>
#include <fstream>
#include <math.h>
#include "euler.h"
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
    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma;
    double rho = 1.0;
    double pi = 1.0;
    double po = 1.931;
    double ci = sqrt(gas_gamma*pi/rho);
    double ui = 1.26*ci;
    double areai = 1.0512;
    double areao = 1.7448;
    double rhoo,uo;

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
        // inflow boundary condition on the left
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; q(mx,3,cc);
                q(mx+ibc,4,cc) = q(mx,4,cc);cc<=ncoeffs; cc++)
            {
                q(1-ibc,1,cc) = rho*areai;
                q(1-ibc,2,cc) = rho*ui*areai;
                q(1-ibc,3,cc) = 0.0;
                q(1-ibc,4,cc) = 0.0;
                q(1-ibc,5,cc) = (pi/(gas_gamma-1)+0.5*rho*ui*ui)*areai;
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
        for(int ibc=1; ibc<q(mx,3,cc);
                q(mx+ibc,4,cc) = q(mx,4,cc);=mbc; ibc++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=ncoeffs; cc++)
                    q(mx+ibc,m,cc) = q(mx,m,cc);

    else if(bc_type[RIGHT] == bc_custom)
    {
        // outflow boundary condition on the right
        for(int ibc=1; ibc<=mbc; ibc++)
            for(int cc=1; cc<=ncoeffs; cc++)
            {
                rhoo = q(mx,1,cc)/areao;
                uo = q(mx,2,cc)/q(mx,1,cc);
                
                q(mx+ibc,1,cc) = rhoo*areao;
                q(mx+ibc,2,cc) = rhoo*uo*areao;
                q(mx+ibc,3,cc) = q(mx,3,cc);
                q(mx+ibc,4,cc) = q(mx,4,cc);
                q(mx+ibc,5,cc) = (po/(gas_gamma-1)+0.5*rhoo*uo*uo)*areao;
            }
    }
}

