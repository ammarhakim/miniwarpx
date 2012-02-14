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
  xloc [in]   - Coordinates in physical space at which kappa is requested
  kappa [out] - Capacity function

 */
void
setkappa(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& kappa)
{
    if(rd.has_kappa)
        std::cout << "** MINIWARPX: setkapp is not provided or is not linked properly" << std::endl;

    int mx = rd.mx;
    int mbc = rd.mbc;
    for(int i=1-mbc; i<=mx+mbc; i++)
        kappa(i) = 1.0; // this means there is no capacity function
}
