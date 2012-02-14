#include "miniwarpx.h"

// flux function for Burger's equation
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    for(int i=1-mbc; i<=mx+mbc; i++)
        fx(i,1) = 0.5*q(i,1)*q(i,1);
}
