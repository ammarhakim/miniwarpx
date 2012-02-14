#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, n, E;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = tfv->mi;
    double me = tfv->me;

    double ratio_me_mi = me/mi;
    double ti = 1.e1;
    double te = 100.*ti;
    double xc = 50.0;
    double omega = 10.;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
//        E = pow(2.*omega,0.5)/cosh(pow(omega,0.5)*(xcell-xc));
        E = 1./cosh(xcell-xc);
        n = -E*E+1;

        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = ratio_me_mi*n;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = 0.0;
        q(i,5)  = n*te/gas_gamma1;
        q(i,6)  = n;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        q(i,10) = n*ti/gas_gamma1;
        q(i,11) = E;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = 0.0;
        q(i,16) = 0.0;
    }
}
