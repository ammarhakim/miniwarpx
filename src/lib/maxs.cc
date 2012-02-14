#include "miniwarpx.h"
#include <iostream>

// dummy function: does nothing. The user MUST override this if the
// simulation is to work properly

/*
  Computes the maxmimum wave speed for equation system

  Parameters
  ----------

  rd [in]  - Input data for simulation
  s  [out] - Maximum wave speed. s(i) is the wave speed in cell i
  q  [in]  - Conserved variable at which speed is to be computed

  Notes
  -----

  This routine needs to be provided only for the MacCormick
  scheme. The maximum wave speed need not be exact: any resonable
  estimate will do.

*/

void 
maxs(const Run_Data& rd, FArray<double>& s, FArray<double>& q)
{
    std::cout << "** MINIWARPX: maxs is not provided or is not linked properly" << std::endl;
}
