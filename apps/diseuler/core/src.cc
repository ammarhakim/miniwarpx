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
    int mx = rd.mx;

    double rhou,rhov,rhow;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double qbym;
    qbym = ev->qbym;
    double gas_gamma = ev->gas_gamma;

    // compute magnetic field
    //magfield(rd,rd.xcoords,ev->bf);

    for (int i=1; i<=mx; i++)
    {
        // compute primitive variables
        rhou   = q(i,2);
        rhov   = q(i,3);
        rhow   = q(i,4);
        
        // compute source terms
        sr(i,1) = 0.0;
        sr(i,2) = qbym*(rhov*ev->bf(i,3)-rhow*ev->bf(i,2));
        sr(i,3) = qbym*(rhow*ev->bf(i,1)-rhou*ev->bf(i,3));
        sr(i,4) = qbym*(rhou*ev->bf(i,2)-rhov*ev->bf(i,1));
        sr(i,5) = 0.0;
        
    }
    
    if(ev->is_radial)
    {
        double r,rho,u,v,w,p;

        for (int i=1; i<=mx; i++)
        {
            r = rd.xcoords(i); // radial coordinate

            // compute primitive variables
            rho  = q(i,1);
            u = q(i,2)/rho;
            v = q(i,3)/rho;
            w = q(i,4)/rho;
            p = (gas_gamma-1)*(q(i,5) - 0.5*rho*(u*u + v*v + w*w));

            // compute source terms
            sr(i,1) += -rho*u/r;
            sr(i,2) += -rho*u*u/r + rho*v*v/r;
            sr(i,3) += -2.*rho*u*v/r;
            sr(i,4) += -rho*u*w/r;
            sr(i,5) += -u*(q(i,5)+p)/r;
        }
    }
}
