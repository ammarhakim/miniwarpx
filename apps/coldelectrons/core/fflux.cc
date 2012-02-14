#include "miniwarpx.h"
#include "cold.h"

// flux function for pressureless gas equation
void 
fflux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    
    double rho, u, v, w;

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
    }
}
