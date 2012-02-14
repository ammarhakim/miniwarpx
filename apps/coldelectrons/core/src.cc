#include "miniwarpx.h"
#include "cold.h"

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
void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double rhoe,rhoeu,rhoev,rhoew;
    double ex,ey,ez,bx,by,bz;

    Cold_Vars *cv = (Cold_Vars*) rd.mvar;

    //  define charge to mass ratio of electrons and ions
    double qbyme = cv->q/cv->m;
    double epsilon0 = cv->epsilon0;

    for (int i=1-mbc; i<=mx+mbc; ++i)
    {
        // compute primitive variables

        // electron fluid
        rhoe   = q(i,1);
        rhoeu  = q(i,2);
        rhoev  = q(i,3);
        rhoew  = q(i,4);
        // electromagnetic field
        ex     = q(i,5);
        ey     = q(i,6);
        ez     = q(i,7);
        bx     = q(i,8);
        by     = q(i,9);
        bz     = q(i,10);

        // electron fluid source terms
        sr(i,1) = 0.0;
        sr(i,2) = qbyme*(rhoe*ex+rhoev*bz-rhoew*by);
        sr(i,3) = qbyme*(rhoe*ey+rhoew*bx-rhoeu*bz);
        sr(i,4) = qbyme*(rhoe*ez+rhoeu*by-rhoev*bx);
 
        // electric field source terms
        sr(i,5) = -qbyme*rhoeu/epsilon0;
        sr(i,6) = -qbyme*rhoev/epsilon0;
        sr(i,7) = -qbyme*rhoew/epsilon0;

        // magnetic field source terms
        sr(i,8) = 0.0;
        sr(i,9) = 0.0;
        sr(i,10) = 0.0;
    }    
}
