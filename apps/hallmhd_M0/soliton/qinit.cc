#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "hallmhd.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, n, b;

    Hallmhd_Vars *hv = (Hallmhd_Vars*) rd.mvar;

    double gas_gamma  = hv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    double ti = 1.e-2;
    double te = 100.*ti;
    double xc = 6.0;

    hv->ur    = 1.; //sqrt(ti);
    hv->betar = 1.;
    hv->rli   = 1.;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        n = 1.0+exp(-0.5*32*fabs(xcell-xc));
//      b = 1.0+exp(-0.5*pow(xcell-xc,2));
        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = n;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = 0.0;
        q(i,5) = n*ti/gas_gamma1;
        q(i,6)  = n;
        q(i,7)  = n*te/gas_gamma1;
        q(i,8) = 0.0;
        q(i,9) = 0.0;
        q(i,10) = 0.0;
        q(i,11) = 0.0;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
    }
}
