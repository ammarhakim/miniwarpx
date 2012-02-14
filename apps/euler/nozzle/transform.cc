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

    double x0 = 4.8;       // point of clustering
    double deg = 0;        // degree of clustering
    double num1,num2,R,i0; // parameters to be calculated for clustering
    double a,b,c,d;        // clustering parameters
    double xl = rd.xlower; // domain in physical space is [xl,xu]
    double xu = rd.xupper; // domain in physical space is [xl,xu]

    a    = pow(10.0,-deg);
    num1 = (xu-x0)/a;
    num2 = (xl-x0)/a;
    R    = log(num1+sqrt(pow(num1,2)+1.0))/ log(num2+sqrt(pow(num2,2)+1.0));
    i0   = (mx-R)/(1.0-R);
    b    = log(num2+sqrt(pow(num2,2)+1.0))/(1-i0);
    c    = -b*i0;
    d    = x0;

    // following grid stretching transformation clusters grid points
    // around x=x0

    // loop over all interior points and apply transform
    for(int i=1; i<=mx; i++)
    {
        xp(i) = a*sinh(b*i+c)+d;
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
