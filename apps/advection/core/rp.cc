#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"

void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    double u = rd.rpar[0]; // advection speed
    int mx = rd.mx;
    int mbc = rd.mbc;

    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        wave(i,1,1) = df(i,1);
        s(i,1) = u;
        if(s(i,1)<0)
            amdq(i,1) = wave(i,1,1);
        else
            apdq(i,1) = wave(i,1,1);
    }

}
