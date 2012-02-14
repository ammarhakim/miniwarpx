#include "miniwarpx.h"

// dummy function: does nothing. The user MUST provide this if
// rd.has_kappa flag is defined

/*
  Sets the capacity function for capacity-form differencing. WarpX is
  capable of solving the system of equations

  kappa*dq/dt + div(f) = s

  Here `kappa` is the capacity function and usually shows up for
  problems on mapped grids.

  Parameters
  ----------

  rd [in]     - Input data for simulation
  kappa [out] - Capacity function

 */
void
setkappa(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& kappa)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    double xlower = rd.xlower;
    double dx = rd.dx;

    FArray<double> xc(Range(1-mbc,mx+mbc),0.0);
    FArray<double> xp(Range(1-mbc,mx+mbc),0.0);

    // set xc to left edge coordinates in computational domain
    for(int i=1-mbc; i<=mx+mbc; i++)
        xc(i) = xlower + (i-1.)*dx; 

    // compute left edge coordinates in physical domain
    grid_transform(rd,xc,xp);

    // kappa is ratio of grid spacing in physical domain to that in
    // computational domain

    // interior cells
    for(int i=1; i<=mx; i++)
        kappa(i) = (xp(i+1)-xp(i))/dx;

    // boundary cells
    for(int i=1; i<=mbc; i++)
    {
        // left edge
        kappa(1-i)  = kappa(i);
        // right edge
        kappa(i+mx) = kappa(mx-i+1);
    }
}
