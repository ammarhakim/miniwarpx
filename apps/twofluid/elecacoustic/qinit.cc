#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell,n,ue,u0,Pe,Ex,E0,rhoe,rho0;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = tfv->mi;
    double me = tfv->me;

    double ratio_me_mi = me/mi;
    double te = 1e-5;
    double ti = 100.*te;
    double xstep = 0.5;
    double P0 = 1.0;
    double L = 1.0;
    double pi = 3.14159;

    n = 1.0;
    u0 = 1.e-8;
    P0 = 1.0;
    E0 = 1.0;
    rhoe = me*n;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        for(int j=1; j<=9; j++)
        {
            // Fourier series for a square wave
            ue = u0*4/pi*1/j*sin(j*pi*xcell/L);
            Pe = P0*4/pi*1/j*sin(j*pi*xcell/L);
            Ex = E0*4/pi*1/j*sin(j*pi*xcell/L);
            //rhoe = rho0*4/pi*1/j*sin(j*pi*xcell/L);
        }

        // IC for electron fluid
        //q(i,1)  = ratio_me_mi*n;
        q(i,1)  = rhoe;
        q(i,2)  = rhoe*ue;
        q(i,3)  = 0.0;
        q(i,4)  = 0.0;
        //q(i,5)  = gas_gamma1*n*te;
        q(i,5)  = Pe/gas_gamma1;
        // IC for ion fluid
        q(i,6)  = n;
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
