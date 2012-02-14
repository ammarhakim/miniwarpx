#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, n;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;

    double gas_gamma  = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    double t = 1.e-2;
    double xc = 0.5;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        n = 1.0+exp(-0.5*32*fabs(xcell-xc));
        // initial conditions for Euler ion fluid
        q(i,1)  = n;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = 0.0;
        q(i,5)  = n*t/gas_gamma1;
    }
}
