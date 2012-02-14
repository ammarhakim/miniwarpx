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

#include "hallmhd.h"

void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    double dx = rd.dx;

    double ne;
    double jx,jy,jz,B;
    FArray<double> gradpex(Range(1-mbc,mx+mbc),0.0);
    double ni,ui,vi,wi;
    double ex,ey,ez,bx,by,bz;

    Hallmhd_Vars *hv = (Hallmhd_Vars*) rd.mvar;

    //  define charge to mass ratio of electrons and ions
    double gas_gamma = hv->gas_gamma;
    double gas_gamma1= gas_gamma-1.;
    double c0 = hv->c0;
    double ur = hv->ur;
    double betar = hv->betar;
    double rli = hv->rli;

    for (int i=1; i<=mx; i++)
    {
        // compute grad Pex (perp to B) using central differencing to get 
        // the cell center values
        if (i>1 && i<mx) gradpex(i) = gas_gamma1*(q(i+1,7)-q(i-1,7))/(2.*dx);
        else if (i==1) gradpex(i) = gas_gamma1*(q(i+2,7)-q(i,7))/(2.*dx);
        else if (i==mx) gradpex(i) = gas_gamma1*(q(i,7)-q(i-2,7))/(2.*dx);
    }

    for (int i=1; i<=mx; i++)
    {
        // compute primitive variables

        // ion fluid
        ni   = q(i,1);
        ui   = q(i,2)/ni;
        vi   = q(i,3)/ni;
        wi   = q(i,4)/ni;
        // electron fluid
        ne   = q(i,6);
        // electromagnetic field
        ex   = q(i,8);
        ey   = q(i,9);
        ez   = q(i,10);
        bx   = q(i,11);
        by   = q(i,12);
        bz   = q(i,13);

        // compute B^2 and jx (jx in perp direction to B)
        B    = sqrt(bx*bx + by*by + bz*bz);
        if (B != 0.) 
        {
            jx   = -(rli*gradpex(i)/B + ne*ex/B - ni*ui);
            jy   =                    -(ne*ey/B - ni*vi);
            jz   =                    -(ne*ez/B - ni*wi);
        }
        else
        {
            jx = ni*ui;
            jy = ni*vi;
            jz = ni*wi;
        }

        // ion fluid source terms
        sr(i,1)  = 0.0;
        sr(i,2)  = ni/rli*(ex+vi*bz-wi*by);
        sr(i,3)  = ni/rli*(ey+wi*bx-ui*bz);
        sr(i,4)  = ni/rli*(ez+ui*by-vi*bx);
        sr(i,5)  = ni*(ui*ex+vi*ey+wi*ez);

        // electron fluid source terms
        sr(i,6) = 0.0;
        sr(i,7) = ex*(jx-ni*ui)+ey*(jy-ni*vi)+ez*(jz-ni*wi);

        // electric field source terms
        sr(i,8) = -betar*c0*c0/(rli*ur*ur)*jx;
        sr(i,9) = -betar*c0*c0/(rli*ur*ur)*jy;
        sr(i,10) = -betar*c0*c0/(rli*ur*ur)*jz;


//        if (B != 0.) 
//        {
            // electric field source terms
//            sr(i,8) = -betar*c0*c0/(rli*B*ur*ur)*(ni*(vi*bz-wi*by)+rli*gradpex(i)+ne*ex);
//            sr(i,9) = -betar*c0*c0/(rli*B*ur*ur)* ni*(wi*bx-ui*bz);
//            sr(i,10) = -betar*c0*c0/(rli*B*ur*ur)* ni*(ui*by-vi*bx);
//        }
//        else
//        {
            // electric field source terms
//            sr(i,8) = 0.;
//            sr(i,9) = 0.;
//            sr(i,10) = 0.;
//        }

        // magnetic field source terms
        sr(i,11) = 0.0;
        sr(i,12) = 0.0;
        sr(i,13) = 0.0;

    }    
}
