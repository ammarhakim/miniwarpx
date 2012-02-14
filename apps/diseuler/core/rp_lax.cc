#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

extern void maxs(const Run_Data& rd, FArray<double>& s, FArray<double>& q);
extern void flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q);

void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    // Reimann solver using Rusonov/Local-Lax fluxes.

    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;
    int mwave = rd.mwave;
    double smax;

    // compute maximum wave speeds
    FArray<double> sl( Range(1-mbc,mx+mbc), 0.0);
    FArray<double> sr( Range(1-mbc,mx+mbc), 0.0);

    maxs(rd,sl,ql); // left  edge speeds
    maxs(rd,sr,qr); // right edge speeds

    // compute waves
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        // 3 waves
        for(int m=1; m<=meqn; m++)
            for(int j=1; j<=mwave; j++)
                wave(i,m,j) = ql(i,m)-qr(i-1,m);
        s(i,1) = dmax(fabs(sl(i)), fabs(sr(i-1)));

        // compute fluctuations
        for(int m=1; m<=meqn; m++)
        {
            apdq(i,m) = s(i,1)*ql(i,m);
            amdq(i,m) = s(i,1)*qr(i-1,m);
        }
    }
}
