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
    // this is a dummy function: it simply sets xp to xc
    int mx = rd.mx;
    int mbc = rd.mbc;

    for(int i=1-mbc; i<=mx+mbc; i++) xp(i) = xc(i);
}
