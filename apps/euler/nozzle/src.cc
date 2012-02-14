#include "miniwarpx.h"
#include <iostream>
#include <math.h>

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
    int mx = rd.mx;
    double xcell, p, area, dadx;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma;

    for (int i=1; i<=mx; i++)
    {
        xcell = rd.xcoords(i);

        // compute pressure and dA/dx to be used in the source term
        area = 1.398+0.347*tanh(0.8*xcell-4.); 
        p    = (gas_gamma-1) * (q(i,5) - 0.5*q(i,2)*q(i,2)/q(i,1))/area;
        dadx = 0.347*0.8 / pow(cosh(0.8*xcell-4.),2);

        // compute source terms
        sr(i,1) = 0.0;
        sr(i,2) = p*dadx;
        sr(i,3) = 0.0;
        sr(i,4) = 0.0;
        sr(i,5) = 0.0;
    }    
}
