#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell,u0,u1,P0,rho0,Bz;
    double nmode,kn,wn,ue,ve,Pe,rhoe;  // Fourier transform variables
    double cs,wc;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;

    double gas_gamma  = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double qbym       = ev->qbym;

    double pi = 3.14159;

    u0 = 1.e-8;
    P0 = 1.0;
    rho0 = 1.0;

    Bz = 1.0;
    cs = sqrt(gas_gamma*P0/rho0);
    wc = qbym*Bz;

    nmode = 9;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        ue  =0;
        rhoe=0;
        ve  =0;
        Pe  =0;
        for(int j=0; j<=nmode; j++)
        {
            // Fourier series for a square wave
            kn   = 2*pi*(2*j+1);
            wn   = sqrt(kn*kn*cs*cs + wc*wc);
            u1   = -u0/(2*j+1)*sin(kn*xcell);
            ue   = ue  + u1;

            rhoe = rhoe-kn/wn*rho0 * u1;
            ve   = ve  -wc/wn*rho0 * u0/(2*j+1)*cos(kn*xcell);
            Pe   = Pe-gas_gamma*P0*kn/wn * u1;
        }
        rhoe = rho0+rhoe;
        Pe   = P0+Pe;

        // Initial conditions
        q(i,1)  = rhoe;
        q(i,2)  = rhoe*ue;
        q(i,3)  = rhoe*ve;
        q(i,4)  = 0.0;
        q(i,5)  = Pe/gas_gamma1+0.5*rho0*(ue*ue+ve*ve);
   }
    for(int i=1-rd.mbc; i<=rd.mx+rd.mbc; i++)
    {
        xcell = xloc(i);

        ev->bf(i,1) = 0.0;
        ev->bf(i,2) = 0.0;
        ev->bf(i,3) = Bz;
    }

}
