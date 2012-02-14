#include <stdio.h>
#include <iostream>
#include <math.h>
#include "cold.h"
#include "miniwarpx.h"

/*
  Reimann solver for equation system.

  Parameters
  ----------
  
  rd [in] - Input data for simulation
  ql [in] - Conserved variables evaluated at left cell edge. Thus ql(i,*,*) is the 
            conserved variable at the left edge of cell i
  qr [in] - Conserved variables evaluated at right cell edge. Thus qr(i,*,*) is the 
            conserved variable at the right edge of cell i
  df [in] - Vector to split at each cell edge. This is usuall the jump in q at each
            edge or the jump in flux at each edge. However, for the RKDG algorithm it
            could be some other vector too. The user should not make any assumptions
            about it in this routine
  wave [out] - wave(i,m,mw) is the mth component of wave mw at edge i.
  s    [out] - s(i,mw) is the speed of wave mw
  amdq [out] - amdq(i,m) is the mth component of the left-going fluctuation at edge i
  apdq [out] - apdq(i,m) is the mth component of the right-going fluctuation at edge i


  Notes
  -----

  For definition of waves and fluctuations see LeVeque's book or WarpX
  notes. For a cell i its left edge is labeled i and its right edge is
  labelled i+1. Thus there is one extra cell edge than cells. For
  standard Roe fluxes the routine eval_fluctuations() can be
  used. This automatically takes care of the slight difference in the
  definition of the fluctuations depending on if the jump in q
  (q-waves) or the jump in flux (f-wave) approach is being used.

*/
void 
rp_electrons(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
             FArray<double>& wave, FArray<double>& s, 
             FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;

    FArray<double> f( Range(1-mbc,mx+mbc, 1,meqn), 0.0 );
    double rl, rr, ul, ur, uav;

    for(int i=2-mbc; i<=mx+mbc; i++)
    {

        ul = qr(i-1,2)/qr(i-1,1); // speed in left cell
        ur = ql(i,2)/ql(i,1); // speed in right cell

        if((ul < 0) && (0 < ur))
        { // vacuum intermediate state will be formed
            fflux(rd,f,qr);
            for(int m=1; m<=meqn; ++m)
                wave(i,m,1) = -f(i-1,m);
            s(i,1) = ul;

            //fflux(rd,f,ql);
            for(int m=1; m<=meqn; ++m)
                wave(i,m,2) = f(i,m);
            s(i,2) = ur;
        }
        else
        { // no vacumm state
            rl = qr(i-1,1);
            rr = ql(i,1);
            // compute Roe avaeraged speed
            uav = (sqrt(rl)*ul + sqrt(rr)*ur)/(sqrt(rl)+sqrt(rr));
            
            if(uav<0)
            {
                for(int m=1; m<=meqn; ++m)
                    wave(i,m,1) = df(i,m);

                for(int m=1; m<=meqn; ++m)
                    wave(i,m,2) = 0.0;
            }
            else
            {
                for(int m=1; m<=meqn; ++m)
                    wave(i,m,1) = 0;

                for(int m=1; m<=meqn; ++m)
                    wave(i,m,2) = df(i,m);
            }
            s(i,1) = uav;
            s(i,2) = uav;
        }
    }
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);
}
