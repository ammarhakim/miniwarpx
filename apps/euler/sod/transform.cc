#include <math.h>
#include "miniwarpx.h"

/*
  Transform the coordinates in the computational domain to physical
  domain.

  Parameters
  ----------

  rd [in]  - Run_Data object
  xc [in]  - Coordinates in computational domain 
  xp [out] - Coordinates in physical domain

  The range of xc and xp is Range(1-mx,mx+mbc)

 */
void
grid_transform(const Run_Data& rd, const FArray<double>& xc, FArray<double>& xp)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double x0 = 0.5;
    double r = 10.0;
    double B,eta;
    double h = rd.xupper; // domain in physical space is [0,h]

    B = log( (1+(exp(r)-1)*x0/h)/(1-(1-exp(-r))*x0/h) )/(2.*r);

    // following grid stretching transformation clusters grid points
    // around x=x0

    // loop over all interior points and apply transform
    for(int i=1; i<=mx; i++)
    {
        eta = xc(i);
        xp(i) = x0*(1 + sinh(r*(eta-B))/sinh(r*B));
    }

    // set ghost cell coordinates by reflection
    for(int i=1; i<=mbc; i++)
    {
        // left edge
        xp(1-i) =  rd.xlower - xp(i);
        // right edge
        xp(i+mx) = rd.xupper + (rd.xupper-xp(mx-i+1));
    }
}
