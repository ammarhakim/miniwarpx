#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double dme = tfv->me;
    double dmi = tfv->mi;
    double dqe = tfv->qe;
    double re = dqe/dme;

    double a,dJ0,dalpha,p,p0,er,djze;
    double pe,pe8,bi,pressure_e,pressure_i,rho_e,rho_i;

    a = 0.125;
    dJ0 = 0.1;
    dalpha = 0.1;
    double r;

    for(int i=1; i<=rd.mx; i++)
    {
        r = rd.xcoords(i); // radial coordinate

        if(r<a) 
        {
            pe = pow(dJ0,2)*(1./4.*pow(r,2) - 12.*pow(r,4) + 4./3.*128.*pow(r,6));
            pe8 = pow(dJ0,2)*(1./4.*pow((1./8.),2) - 12.*pow((1./8.),4)
                              + 4./3.*128.*pow((1./8.),6));
            p0 = pe8/(1-dalpha);
            p = p0 - pe;
            bi = dJ0*(0.5*r-16.*pow(r,3));

            djze = dJ0*(1-64.0*pow(r,2));
            pressure_e = 0.5*p;
            pressure_i = 0.5*p;
            rho_e = dme*p/p0;
            rho_i = dmi*p/p0;
        }
        else
        {
            pe = pow(dJ0,2)*(1./4.*pow((1./8.),2) - 12.*pow((1./8.),4)
                             + 4./3.*128.*pow((1./8.),6));
            p0 = pe/(1-dalpha);
            p = p0-pe;
            bi = dJ0*(0.5*(1./8.)-16.*pow((1./8.),3))*(1./8.)*(1./r);
                
            pressure_e = 0.5*p;
            pressure_i = 0.5*p;
            rho_e = dme*p/p0;
            rho_i = dmi*p/p0;
            djze = 0.0;
        }

        // set electron variables
        q(i,1) = rho_e;
        q(i,2) = 0;
        q(i,3) = 0;
        q(i,4) = djze/re;
        er = pressure_e/gas_gamma1 + 0.5*pow(q(i,4),2)/rho_e;
        q(i,5) = er;
        // set ion variables
        q(i,6) = rho_i;
        q(i,7) = 0;
        q(i,8) = 0;
        q(i,9) = 0;
        er = pressure_i/gas_gamma1;
        q(i,10) = er;
        // EM fields
        q(i,11) = 0;
        q(i,12) = 0;
        q(i,13) = 0;
        q(i,14) = 0;
        q(i,15) = bi;
        q(i,16) = 0;
    }

}
