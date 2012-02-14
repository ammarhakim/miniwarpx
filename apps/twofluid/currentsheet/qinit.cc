#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell,p,vi,ve,mey,miy,ne,ni,pe,pi,Bz;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double epsilon0 = tfv->epsilon0;
    double mi = tfv->mi;
    double me = tfv->me;
    double qi = tfv->qi;
    double qe = tfv->qe;
    double c0 = tfv->c0;
    double Lx = 2.0;
    double Ti = 0.1;
    double Te = 0.1;
    double B0 = pow(2*(Ti+Te),0.5);

    double ratio_me_mi = me/mi;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        ne = pow(1/(cosh(xcell/Lx)),2.0);
        ni = ne;
        vi = -2*c0*Ti/(qi*B0*Lx);
        ve = -vi*(Te/Ti);
        mey = ne*me*ve;
        miy = ni*mi*vi;
        pe = Te*pow(1/(cosh(xcell/Lx)),2.0);
        pi = Ti*pow(1/(cosh(xcell/Lx)),2.0);
        Bz = B0*tanh(xcell/Lx);

        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = ne*me;
        q(i,2)  = 0.0;
        q(i,3)  = mey;
        q(i,4)  = 0.0;
        q(i,5)  = pe/gas_gamma1 + 0.5*pow(mey,2)/(ne*me);
        q(i,6)  = ni*mi;
        q(i,7)  = 0.0;
        q(i,8)  = miy;
        q(i,9)  = 0.0;
        q(i,10) = pi/gas_gamma1 + 0.5*pow(miy,2)/(ni*mi);
        q(i,11) = 0.0;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = 0.0;
        q(i,16) = Bz;
    }
}
