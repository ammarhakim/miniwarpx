#include "miniwarpx.h"
#include <iostream>

/*
  Computes the maxmimum wave speed tenmoment equations

  Parameters
  ----------

  rd [in]  - Input data for simulation
  s  [out] - Maximum wave speed. s(i) is the wave speed in cell i
  q  [in]  - Conserved variable at which speed is to be computed

*/

void 
maxs(const Run_Data& rd, FArray<double>& s, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    


    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
    }
}
