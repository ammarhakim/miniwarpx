#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

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

  Euler Note
  ----------

  This Reimann solver is for the complete Euler equation system (5
  equations). There are 5 eigenvalues: u-c,u,u,u,u+c, three of which
  are lumped into a single wave. Thus mwave=3 for this routine.

*/
void 
rp_euler(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
         FArray<double>& wave, FArray<double>& s, 
         FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double rhsqrtl,rhsqrtr,pl,pr,rhsq2;
    double a1,a2,a3,a4,a5;

    FArray<double> delta(Range(1,5),0.0);

    // get hold of Two-Fluid specific data
    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;
    // gas constant
    double gas_gamma = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    // compute Roe-averaged quantities
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        if ((qr(i-1,1)<0) || (ql(i,1)<0))
        {
            std::cout << "*** Negative density in routine rp in file rp3eu.cc" << std::endl;
            exit(1); // abort execution
        }

        rhsqrtl = sqrt(qr(i-1,1));
        rhsqrtr = sqrt(ql(i,1));

        // left edge pressure
        pl = gas_gamma1*(qr(i-1,5) 
                         - 0.5*(pow(qr(i-1,2),2) +
                                pow(qr(i-1,3),2) +
                                pow(qr(i-1,4),2))/qr(i-1,1));
        // right edge pressure
        pr = gas_gamma1*(ql(i,5) 
                         - 0.5*(pow(ql(i,2),2) +
                                pow(ql(i,3),2) +
                                pow(ql(i,4),2))/ql(i,1));

        if ((pl<0) || (pr<0))
        {
            std::cout << "*** Negative pressure in routine rp in file rp3eu.cc" << std::endl;
            exit(1); // abort execution
        }

        
        rhsq2 = rhsqrtl + rhsqrtr;

        // Roe-averaged velocity components
        tfv->u(i) = (qr(i-1,2)/rhsqrtl + ql(i,2)/rhsqrtr) / rhsq2;
        tfv->v(i) = (qr(i-1,3)/rhsqrtl + ql(i,3)/rhsqrtr) / rhsq2;
        tfv->w(i) = (qr(i-1,4)/rhsqrtl + ql(i,4)/rhsqrtr) / rhsq2;

        tfv->u2v2w2(i) = pow(tfv->u(i),2) + pow(tfv->v(i),2) + pow(tfv->w(i),2);

        // Roe-averaged enthalpy
        tfv->enth(i) = (((qr(i-1,5)+pl)/rhsqrtl
                        + (ql(i,5)+pr)/rhsqrtr)) / rhsq2;

        // speed of sound
        a2 = gas_gamma1*(tfv->enth(i) - .5*tfv->u2v2w2(i));
        tfv->a(i) = sqrt(a2);

        if(a2<0)
        {
            std::cout << "*** Negative sound-speed in routine rp in file rp_euler.cc" << std::endl;
            exit(1); // abort execution
        }

        tfv->g1a2(i) = gas_gamma1 / a2;
        tfv->euv(i)  = tfv->enth(i) - tfv->u2v2w2(i);
    }

    // compute waves
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        // compute coefficients of the 5 eigenvectors
        delta(1) = df(i,1);
        delta(2) = df(i,2);
        delta(3) = df(i,3);
        delta(4) = df(i,4);
        delta(5) = df(i,5);

        a4 = tfv->g1a2(i) * (tfv->euv(i)*delta(1)
                            + tfv->u(i)*delta(2) 
                            + tfv->v(i)*delta(3) 
                            + tfv->w(i)*delta(4)
                            - delta(5));
        a2 = delta(3) - tfv->v(i)*delta(1);
        a3 = delta(4) - tfv->w(i)*delta(1);
        a5 = (delta(2) + (tfv->a(i)-tfv->u(i))*delta(1) - tfv->a(i)*a4) / (2.0*tfv->a(i));
        a1 = delta(1) - a4 - a5;

        // Wave 1: eigenvalue is u-c
        wave(i,1,1)  = a1;
        wave(i,2,1) = a1*(tfv->u(i)-tfv->a(i));
        wave(i,3,1) = a1*tfv->v(i);
        wave(i,4,1) = a1*tfv->w(i);
        wave(i,5,1)  = a1*(tfv->enth(i) - tfv->u(i)*tfv->a(i));
        s(i,1) = tfv->u(i)-tfv->a(i);

        // Wave 2: the 3 eigenvectors corresponding to the repeated
        // eigenvalue u,u,u are lumped together into a single wave
        wave(i,1,2)  = a4;
        wave(i,2,2) = a4*tfv->u(i);
        wave(i,3,2) = a4*tfv->v(i)	 	 + a2;
        wave(i,4,2) = a4*tfv->w(i)	 	 + a3;
        wave(i,5,2)  = a4*0.5*tfv->u2v2w2(i)  + a2*tfv->v(i) + a3*tfv->w(i);
        s(i,2) = tfv->u(i);

        // Wave 3: eigenvalue is u+c
        wave(i,1,3)  = a5;
        wave(i,2,3) = a5*(tfv->u(i)+tfv->a(i));
        wave(i,3,3) = a5*tfv->v(i);
        wave(i,4,3) = a5*tfv->w(i);
        wave(i,5,3)  = a5*(tfv->enth(i)+tfv->u(i)*tfv->a(i));
        s(i,3) = tfv->u(i)+tfv->a(i);
    }

    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);

}
