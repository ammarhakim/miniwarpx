#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "mhd.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double n0, B0, pr;

    n0 = 1.0;
    B0 = 0.1;

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma; // gas constant
    double gas_gamma1 = gas_gamma - 1;

    // Harris current sheet equilibium
    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        q(i,1) = 0.5 + n0/pow(cosh(xcell),2.0);
        q(i,2) = 0.0;
        q(i,3) = 0.0;
        q(i,4) = 0.0;

        q(i,6) = 0.0;
        q(i,7) = 0.0;
        q(i,8) = B0*tanh(xcell);
        
        pr = 3.0 - 0.5*pow(q(i,8),2);
        q(i,5) = pr/gas_gamma1 + 0.5*pow(q(i,8),2.);
    }

}
