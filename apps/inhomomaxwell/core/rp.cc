#include <stdio.h>
#include <iostream>
#include <cmath>
#include "miniwarpx.h"
#include "maxwell.h"

void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    double a1, a2;
    double etaim, etai, cim, ci;
    IMaxwell_Vars *v = (IMaxwell_Vars*) rd.mvar;

    FArray<double> delta(Range(1,2));


    // compute waves
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        delta(1) = df(i,1);
        delta(2) = df(i,2);

        // impedances and speeds
        etaim = std::sqrt(v->mu(i-1)/v->ep(i-1));
        etai = std::sqrt(v->mu(i)/v->ep(i));

        cim = 1.0/std::sqrt(v->mu(i-1)*v->ep(i-1));
        ci = 1.0/std::sqrt(v->mu(i)*v->ep(i));
        

        // compute coefficients
        a1 = (-delta(1) + etai*delta(2))/(etaim+etai);
        a2 = (delta(1) + etaim*delta(2))/(etaim+etai);

        // wave 1 with speed -cim
        wave(i,1,1) = -a1*etaim;
        wave(i,2,1) = a1;
        s(i,1) = -cim;

        // wave 2 with speed ci
        wave(i,1,2) = a2*etai;
        wave(i,2,2) = a2;
        s(i,2) = ci;
    }
    
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);
}
