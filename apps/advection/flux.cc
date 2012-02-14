#include "miniwarpx.h"

// flux function for linear-advection equation with constant advection
// velocity
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    double u = rd.rpar[0];
    int mx = rd.mx;
    int mbc = rd.mbc;    

    for(int i=1-mbc; i<=mx+mbc; i++)
        fx(i,1) = u*q(i,1);
}
