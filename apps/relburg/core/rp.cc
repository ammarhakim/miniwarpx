#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"

void
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    double p, ul, ur;

    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        wave(i,1,1) = df(i,1);
        // compute speeds
        p = ql(i,1);
        ul = p/sqrt(1+p*p);
        p = qr(i-1,1);
        ur = p/sqrt(1+p*p);

        s(i,1) = 0.5*(ul+ur);
    }
    eval_fluctuations(rd,wave,s,amdq,apdq);
}
