#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"
#include "tenmom.h"

void 
rp_maxwell(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
           FArray<double>& wave, FArray<double>& s, 
           FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;

    double a1,a2,a3,a4,a5,a6;

    FArray<double> delta(Range(1,6));

    // speed of light
    Tenmom_Vars *tfv = (Tenmom_Vars*) rd.mvar;
    double c0 = tfv->c0;
    
    // compute waves
    for(int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute coefficients of the 4 eigenvectors
        delta(1) = df(i,1);
        delta(2) = df(i,2);
        delta(3) = df(i,3);
        delta(4) = df(i,4);
        delta(5) = df(i,5);
        delta(6) = df(i,6);

        a1 = 0.5*(delta(3)/c0+delta(5));
        a2 = 0.5*(-delta(2)/c0+delta(6));
        a3 = delta(4);
        a4 = delta(1);
        a5 = 0.5*(-delta(3)/c0+delta(5));
        a6 = 0.5*(delta(2)/c0+delta(6));

        // compute waves

        // wave 1 with speed -c
        wave(i,1,1) = 0;
        wave(i,2,1) = -a2*c0;
        wave(i,3,1) = a1*c0;
        wave(i,4,1) = 0;
        wave(i,5,1) = a1;
        wave(i,6,1) = a2;
        s(i,1) = -c0;

        // wave 2 with speed 0
        wave(i,1,2) = a4;
        wave(i,2,2) = 0;
        wave(i,3,2) = 0;
        wave(i,4,2) = a3;
        wave(i,5,2) = 0;
        wave(i,6,2) = 0;
        s(i,2) = 0;

        // wave 3 with speed +c
        wave(i,1,3) = 0;
        wave(i,2,3) = a6*c0;
        wave(i,3,3) = -a5*c0;
        wave(i,4,3) = 0;
        wave(i,5,3) = a5;
        wave(i,6,3) = a6;
        s(i,3) = c0;
    }
    
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);
}
