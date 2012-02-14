#include "miniwarpx.h"
#include "strainwave.h"
#include <iostream>

/*
  Computes the maxmimum wave speed 

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
    double dsig_dstrain, c0;

    Strainwave_Vars *sv = (Strainwave_Vars*) rd.mvar;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        dsig_dstrain = sv->modul(i) 
            + 2*sv->beta*sv->modul(i)*sv->modul(i)*q(i,1);
        c0 = sqrt(dsig_dstrain/sv->rho(i));
        s(i) = c0; // maximum wave speed 

    }
}
