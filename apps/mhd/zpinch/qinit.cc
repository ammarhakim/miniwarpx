#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "mhd.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double pr, p0, n0, alpha, beta, ex;

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma; // gas constant
    double gas_gamma1 = gas_gamma - 1;

    p0 = 1.0;
    n0 = 1.0;
    alpha = 1.0;
    beta = 0.1;

    // Smooth Z-pinch profile
    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);
        ex = exp(alpha*xcell*xcell);

        q(i,1) = n0/ex + n0*beta;
        q(i,2) = 0.0;
        q(i,3) = 0.0;
        q(i,4) = 0.0;

        q(i,6) = 0.0;
        q(i,7) = sqrt(2*p0/alpha*(1-(1+alpha*xcell*xcell)/ex))/xcell;
        q(i,8) = 0.0;
        
        pr = p0/ex + p0*beta;
        q(i,5) = pr/gas_gamma1 + 0.5*pow(q(i,7),2.);
    }
}
