#include "miniwarpx.h"
#include <iostream>

// dummy function: does nothing. The user MUST override this if the
// equation system has a source terms.

/*
  Computes the source function for equation system

  Parameters
  ----------

  rd [in]  - Input data for simulation
  sr [out] - source array. sr(i,*) is the source computed using conserved variables q(i,*)
  q  [in]  - Conserved variable at which source is to be computed

  Notes
  -----

  In most cases the source function does not explicity depend on time
  or spatial coordinates. However, when this function is called
  rd.tcurrent is set to the current time and rd.xcoords(i) is set to
  the coordinate in cell i at which the source should be
  computed. Thus, using these, one can compute the sources which
  explicity depend on time and spatial coordinates.
*/

#include "euler.h"

void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    for (int i=1; i<=rd.mx; i++)
    {
        // compute source terms
        sr(i,1) = 0.0;
        sr(i,2) = -0.1*q(i,1);
        sr(i,3) = 0.0;
        sr(i,4) = 0.0;
        sr(i,5) = -0.1*q(i,2);
    }
}
