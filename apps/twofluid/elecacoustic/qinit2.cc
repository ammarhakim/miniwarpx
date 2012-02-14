#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell,n,ue,u0,Pe,Ex,E0,rhoe,rho0,rhoi;
    double kn,wn,cs,wc,u1;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = tfv->mi;
    double me = tfv->me;
    double qe = tfv->qe;
    double epsilon0 = tfv->epsilon0;

    double ratio_me_mi = me/mi;
    double te = 1e-5;
    double ti = 100.*te;
    double xstep = 0.5;
    double P0 = 1.0;
    double L = 1.0;
    double pi = 3.14159;
    int nmode = 9;

    n = 1.0;
    u0 = 1.e-1;
    P0 = 1.0;
    E0 = 1.0;
    rho0 = 1.0;
    rhoe = me*n;
    rhoi = mi*n;

    cs = sqrt(gas_gamma*P0/rho0);
    wc = sqrt(n*qe*qe/(epsilon0*me));

    for(int i=1; i<=rd.mx; i++)
    {

        xcell = xloc(i);
        ue  =0;
        rhoe=0;
        Pe  =0;
        Ex  =0;
        for(int j=0; j<=nmode; j++)
        {
            // Fourier series for a square wave
            kn   = 2*pi*(2*j+1);
            wn   = sqrt(kn*kn*cs*cs + wc*wc);
            u1   = -u0/(2*j+1)*sin(kn*xcell);
            ue   = ue  + u1;

            rhoe = rhoe-kn/wn*rho0 * u1;
            Pe   = Pe-gas_gamma*P0*kn/wn * u1;
            Ex   = Ex-gas_gamma*E0*kn/wn * u1;
        }
        rhoe = rho0+rhoe;
        Pe   = P0+Pe;

        // IC for electron fluid
        //q(i,1)  = ratio_me_mi*n;
        q(i,1)  = rhoe;
        q(i,2)  = rhoe*ue;
        q(i,3)  = 0.0;
        q(i,4)  = 0.0;
        //q(i,5)  = gas_gamma1*n*te;
        q(i,5)  = Pe/gas_gamma1 + 0.5*pow(q(i,2),2)/rhoe;
        // IC for ion fluid
        q(i,6)  = rhoi;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        //q(i,10) = gas_gamma1*n*ti;
        q(i,10)  = P0/gas_gamma1;
        // IC for EM terms
        q(i,11) = Ex;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = 0.0;
        q(i,16) = 0.0;
    }
}
