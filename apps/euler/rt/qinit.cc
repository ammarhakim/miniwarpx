#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, rho, pr;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        if(xcell < 0.0)
            rho = 1.0;
        else
            rho = 1.0;
        pr = 2.0 - xcell*rho*0.1;
        q(i,1) = rho;
        q(i,2) = 0.0;
        q(i,3) = 0.0;
        q(i,4) = 0.0;
        q(i,5) = pr/(gas_gamma-1);
    }
}
