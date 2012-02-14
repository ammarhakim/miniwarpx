#include "miniwarpx.h"

// flux function for pressureless gas equation
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    
    double rho, u;

    for(int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        rho = q(i,1);
        u = q(i,2)/rho;
        
        fx(i,1) = rho*u;
        fx(i,2) = rho*u*u;
    }
}
