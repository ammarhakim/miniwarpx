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
    double nr, nl, ur, ul, ns;

    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        nr = ql(i,1);   ur = ql(i,2);
        nl = qr(i-1,1); ul = qr(i-1,2);
        ns = ul*nl/nr; // intermediate density

        // wave 1
        wave(i,1,1) = (nr-ns)*ur;
        wave(i,2,1) = 0.0;
        s(i,1) = ur;

        // wave 2
        wave(i,1,2) = 0.0;
        wave(i,2,2) = df(i,2);
        s(i,2) = 0.5*(ql(i,2) + qr(i-1,2));
    }
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);    
}
