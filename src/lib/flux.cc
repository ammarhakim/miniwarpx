#include "miniwarpx.h"
#include <iostream>

// dummy function: does nothing. The user MUST override this if the
// simulation is to work properly

/*
  Computes the flux function for equation system

  Parameters
  ----------

  rd [in]  - Input data for simulation
  fx [out] - Flux array. f(i,*) is the flux computed using conserved variables q(i,*)
  q  [in]  - Conserved variable at which flux is to be computed

  Notes
  -----

  In most cases the flux function does not explicity depend on time or
  spatial coordinates. However, when this function is called
  rd.tcurrent is set to the current time and rd.xcoords(i) is set to
  the coordinate in cell i at which the flux should be computed. Thus,
  using these, one can compute the fluxes which explicity depend on
  time and spatial coordinates.
*/

void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    std::cout << "** MINIWARPX: flux is not provided or is not linked properly" << std::endl;
}
