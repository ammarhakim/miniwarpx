#include "miniwarpx.h"
#include "cold.h"

// flux function for pressureless gas equation
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    
    double rho, u, v, w;

    Cold_Vars *cv = (Cold_Vars*) rd.mvar;
    // speed of light
    double c0 = cv->c0;

    for(int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        rho = q(i,1);
        u = q(i,2)/rho;
        v = q(i,3)/rho;
        w = q(i,4)/rho;
        
        // Fluid fluxes
        fx(i,1) = rho*u;
        fx(i,2) = rho*u*u;
        fx(i,3) = rho*u*v;
        fx(i,4) = rho*u*w;
        // EM field fluxes
        fx(i,5) = 0.0;
        fx(i,6) = c0*c0*q(i,10);
        fx(i,7) = -c0*c0*q(i,9);
        fx(i,8) = 0.0;
        fx(i,9) = -q(i,7);
        fx(i,10) = q(i,6);
    }
}
