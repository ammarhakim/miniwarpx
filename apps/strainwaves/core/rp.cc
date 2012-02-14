#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"
#include "strainwave.h"

void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double a1,a2,c0,Z,q1avg,q2avg,dsig_dstrain;

    FArray<double> delta(Range(1,2));

    Strainwave_Vars *sv = (Strainwave_Vars*) rd.mvar;
        
    // compute waves
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        q1avg = (ql(i,1)+qr(i-1,1))/2;
        q2avg = (ql(i,2)+qr(i-1,2))/2;

        // for linear problems beta=0, so can use these statements
        // for both linear and nonlinear problems
        dsig_dstrain = sv->modul(i) 
            + 2*sv->beta*sv->modul(i)*sv->modul(i)*q1avg;
        c0 = sqrt(dsig_dstrain/sv->rho(i));
        Z  = sqrt(sv->rho(i)*dsig_dstrain);

        // compute coefficients of the 2 eigenvectors
        delta(1) = df(i,1);
        delta(2) = df(i,2);

        a1 = 0.5*(delta(1)+delta(2)/Z);
        a2 = 0.5*(delta(1)-delta(2)/Z);

        // compute waves

        // wave 1 with speed -c
        wave(i,1,1) = a1;
        wave(i,2,1) = a1*Z;
        s(i,1) = -c0;

        // wave 2 with speed +c
        wave(i,1,2) = a2;
        wave(i,2,2) = -a2*Z;
        s(i,2) = +c0;
    }
    
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);
}
