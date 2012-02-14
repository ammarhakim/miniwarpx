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

#include "twofluid.h"

void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double rhoe,rhoeu,rhoev,rhoew;
    double rhoi,rhoiu,rhoiv,rhoiw;
    double ex,ey,ez,bx,by,bz;
    double ue,ve,we,ui,vi,wi;
    double Pe,Pi;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    //  define charge to mass ratio of electrons and ions
    double qbyme = tfv->qe/tfv->me;
    double qbymi = tfv->qi/tfv->mi;
    double epsilon0 = tfv->epsilon0;
    double gas_gamma = tfv->gas_gamma;
    double c0 = tfv->c0;
    int is_radial = tfv->is_radial;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables

        // electron fluid
        rhoe   = q(i,1);
        rhoeu  = q(i,2);
        rhoev  = q(i,3);
        rhoew  = q(i,4);
        // ion fluid
        rhoi   = q(i,6);
        rhoiu  = q(i,7);
        rhoiv  = q(i,8);
        rhoiw  = q(i,9);
        // electromagnetic field
        ex     = q(i,11);
        ey     = q(i,12);
        ez     = q(i,13);
        bx     = q(i,14);
        by     = q(i,15);
        bz     = q(i,16);

        // electron fluid source terms
        sr(i,1) = 0.0;
        sr(i,2) = qbyme*(rhoe*ex+rhoev*bz-rhoew*by);
        sr(i,3) = qbyme*(rhoe*ey+rhoew*bx-rhoeu*bz);
        sr(i,4) = qbyme*(rhoe*ez+rhoeu*by-rhoev*bx);
        sr(i,5) = qbyme*(rhoeu*ex+rhoev*ey+rhoew*ez);

        // ion fluid source terms
        sr(i,6)  = 0.0;
        sr(i,7)  = qbymi*(rhoi*ex+rhoiv*bz-rhoiw*by);
        sr(i,8)  = qbymi*(rhoi*ey+rhoiw*bx-rhoiu*bz);
        sr(i,9)  = qbymi*(rhoi*ez+rhoiu*by-rhoiv*bx);
        sr(i,10) = qbymi*(rhoiu*ex+rhoiv*ey+rhoiw*ez);

        // electric field source terms
        sr(i,11) = -(qbyme*rhoeu+qbymi*rhoiu)/epsilon0;
        sr(i,12) = -(qbyme*rhoev+qbymi*rhoiv)/epsilon0;
        sr(i,13) = -(qbyme*rhoew+qbymi*rhoiw)/epsilon0;

        // magnetic field source terms
        sr(i,14) = 0.0;
        sr(i,15) = 0.0;
        sr(i,16) = 0.0;

        if (is_radial==1)
        {
            ue = rhoeu/rhoe;
            ve = rhoev/rhoe;
            we = rhoew/rhoe;
            Pe = (gas_gamma-1.0)*(q(i,5)-0.5*rhoe*(ue*ue+ve*ve+we*we));
            ui = rhoiu/rhoi;
            vi = rhoiv/rhoi;
            wi = rhoiw/rhoi;
            Pi = (gas_gamma-1.0)*(q(i,10)-0.5*rhoi*(ui*ui+vi*vi+wi*wi));
            
            // electron fluid source terms
            sr(i,1) += -rhoeu/rd.xcoords(i);
            sr(i,2) += -rhoe*ue*ue/rd.xcoords(i) + rhoe*ve*ve/rd.xcoords(i);
            sr(i,3) += -2.*rhoe*ue*ve/rd.xcoords(i);
            sr(i,4) += -rhoe*ue*we/rd.xcoords(i);
            sr(i,5) += -ue*(q(i,5)+Pe)/rd.xcoords(i);
            
            // ion fluid source terms
            sr(i,6) += -rhoiu/rd.xcoords(i);
            sr(i,7) += -rhoi*ui*ui/rd.xcoords(i) + rhoi*vi*vi/rd.xcoords(i);
            sr(i,8) += -2.*rhoi*ui*vi/rd.xcoords(i);
            sr(i,9) += -rhoi*ui*wi/rd.xcoords(i);
            sr(i,10) += -ui*(q(i,10)+Pi)/rd.xcoords(i);

            // electric field source terms
            sr(i,11) += 0.0;
            sr(i,12) += 0.0;
            sr(i,13) += c0*c0*by/rd.xcoords(i);
            
            // magnetic field source terms
            sr(i,14) += 0.0;
            sr(i,15) += 0.0;
            sr(i,16) += -ey/rd.xcoords(i);
        }
    }    
}
