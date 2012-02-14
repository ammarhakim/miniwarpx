#include "miniwarpx.h"
#include <iostream>
#include "tenmoment.h"

/*

  Computes the source function for the 10 moment equation system

  Parameters
  ----------

  rd [in] - Input data for simulation
  sr [out]- source array. sr(i,*) is the source computed using conserved variables
  q  [in] - Conserved variable at which the source is to be computed.

*/

void

src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;

    double rhoe,mxe,mye,mze,Pxxe,Pxye,Pxze,Pyye,Pyze,Pzze;
    double ue,ve,we;

    double rhoi,mxi,myi,mzi,Pxxi,Pxyi,Pxzi,Pyyi,Pyzi,Pzzi;
    double ui,vi,wi;

    double Ex,Ey,Ez,Bx,By,Bz;

    TenMoment_Vars *tmv = (TenMoment_Vars*) rd.mvar;

    // qtom (charge to mass ratio)
    double qtomi, qtome;
    double mi,ni,qi,me,ne,qe;
    double c0,c02,muo;

    qtomi = tmv -> qtomi;
    mi    = tmv -> mi;
    qtome = tmv -> qtome;
    me    = tmv -> me;
    muo   = tmv -> muo;
    c0    = tmv -> c0;
    
    qi    = qtomi*mi;
    qe    = qtome*me;

    c02   = c0*c0;

    for (int i=1; i<=mx; i++)
    {
        // compute primitive variables
        // ion fluid
        rhoi= q(i,1);
        mxi = q(i,2);
        myi = q(i,3);
        mzi = q(i,4);
        Pxxi= q(i,5);
        Pxyi= q(i,6);
        Pxzi= q(i,7);
        Pyyi= q(i,8);
        Pyzi= q(i,9);
        Pzzi= q(i,10);

        // electron fluid
        rhoe= q(i,11);
        mxe = q(i,12);
        mye = q(i,13);
        mze = q(i,14);
        Pxxe= q(i,15);
        Pxye= q(i,16);
        Pxze= q(i,17);
        Pyye= q(i,18);
        Pyze= q(i,19);
        Pzze= q(i,20);

        // e&m fields
        Ex  = q(i,21);
        Ey  = q(i,22);
        Ez  = q(i,23);
        Bx  = q(i,24);
        By  = q(i,25);
        Bz  = q(i,26);

        // extracting primitives
        ni = rhoi/mi;
        ui = mxi/rhoi;
        vi = myi/rhoi;
        wi = mzi/rhoi;

        ne = rhoe/me;
        ue = mxe/rhoe;
        ve = mye/rhoe;
        we = mze/rhoe;

        
        // compute source terms
        sr(i,1) =   0.0;
        sr(i,2) =   qtomi* rhoi*(Ex + vi*Bz - wi*By);
        sr(i,3) =   qtomi* rhoi*(Ey + wi*Bx - ui*Bz);
        sr(i,4) =   qtomi* rhoi*(Ez + ui*By - vi*Bx);
        sr(i,5) = 2*qtomi*(rhoi* ui*Ex + Bz*Pxyi - By*Pxzi);
        sr(i,6) =   qtomi* rhoi*(ui*Ey + vi*Ex) + qtomi*(Bz*Pyyi - By*Pyzi - Bz*Pxxi + Bx*Pxzi);
        sr(i,7) =   qtomi* rhoi*(ui*Ez + wi*Ex) + qtomi*(Bz*Pyzi + By*Pxxi - By*Pzzi - Bx*Pxyi);
        sr(i,8) = 2*qtomi*(rhoi* vi*Ey + Bx*Pyzi - Bz*Pxyi);
        sr(i,9) =   qtomi* rhoi*(vi*Ez + wi*Ey) + qtomi*(By*Pxyi - Bz*Pxzi + Bx*Pzzi - Bx*Pyyi);
        sr(i,10)= 2*qtomi*(rhoi* wi*Ez + By*Pxzi - Bx*Pyzi);

        sr(i,11) =   0.0;
        sr(i,12) =   qtome* rhoe*(Ex + ve*Bz - we*By);
        sr(i,13) =   qtome* rhoe*(Ey + we*Bx - ue*Bz);
        sr(i,14) =   qtome* rhoe*(Ez + ue*By - ve*Bx);
        sr(i,15) = 2*qtome*(rhoe* ue*Ex + Bz*Pxye - By*Pxze);
        sr(i,16) =   qtome* rhoe*(ue*Ey + ve*Ex) + qtome*(Bz*Pyye - By*Pyze - Bz*Pxxe + Bx*Pxze);
        sr(i,17) =   qtome* rhoe*(ue*Ez + we*Ex) + qtome*(Bz*Pyze + By*Pxxe - By*Pzze - Bx*Pxye);
        sr(i,18) = 2*qtome*(rhoe* ve*Ey + Bx*Pyze - Bz*Pxye);
        sr(i,19) =   qtome* rhoe*(ve*Ez + we*Ey) + qtome*(By*Pxye - Bz*Pxze + Bx*Pzze - Bx*Pyye);
        sr(i,20) = 2*qtome*(rhoe* we*Ez + By*Pxze - Bx*Pyze);

        sr(i,21) = -c02*muo*(ni*qi*ui + ne*qe*ue);
        sr(i,22) = -c02*muo*(ni*qi*vi + ne*qe*ve);
        sr(i,23) = -c02*muo*(ni*qi*wi + ne*qe*we);
        sr(i,24) = 0;
        sr(i,25) = 0;
        sr(i,26) = 0;
    }
}
