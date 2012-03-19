#include "miniwarpx.h"
#include <iostream>

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

#include "maxwell.h"
#include <iostream>

#include <cmath>

void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    Maxwell_Vars *mv = (Maxwell_Vars*) rd.mvar;
    double t = rd.tcurrent;

    double dx = (rd.xupper-rd.xlower)/rd.mx;
    double xLastEdge = rd.xupper - dx;

    double pi = 3.141592654;
    double driveOmega = 23545644591.361;
    double driveF = driveOmega/(2*pi);

    //  speed of light
    double J0 = 1.0;
    double epsilon0 = 8.854187817e-12;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // electron fluid source terms
        sr(i,1) = 0.0;
        if ((rd.xcoords(i) > xLastEdge) && (rd.xcoords(i) < rd.xupper))
        {
          double ramp = std::sin(0.5*pi*std::min(1.0, 0.1*driveF*t));
          sr(i,2) = -J0*ramp*ramp*std::sin(driveOmega*t)/epsilon0;
        }
        else
          sr(i,2) = 0.0;
        sr(i,3) = 0.0;
        sr(i,4) = 0.0;
        sr(i,5) = 0.0;
        sr(i,6)  = 0.0;
    }
}
