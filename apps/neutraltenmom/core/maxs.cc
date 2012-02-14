#include "miniwarpx.h"
#include "tenmoment.h"
#include <iostream>

/*
  Computes the maxmimum wave speed for Two-Fluid equations: this is
  simply speed of light

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

    TenMoment_Vars *tmv = (TenMoment_Vars*) rd.mvar;
    // speed of light
    double c0 = tmv->c0;

    for (int i=1-mbc; i<=mx+mbc; i++)
        s(i) = c0; // maximum wave speed is simply speed of light
}
