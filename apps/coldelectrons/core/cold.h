#ifndef __cold__
#define __cold__

#include "miniwarpx.h"

enum {RHO = 1, RHOU, RHOV, RHOW, EX, EY, EZ, BX, BY, BZ};

struct Cold_Vars
{
    double c0; // speed of light
    double q, m; // electron charge, mass
    double epsilon0; // 
    FArray<double> n0; // background ion-density
};

void 
rp_maxwell(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
           FArray<double>& wave, FArray<double>& s, 
           FArray<double>& amdq, FArray<double>& apdq);
void 
rp_electrons(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
             FArray<double>& wave, FArray<double>& s, 
             FArray<double>& amdq, FArray<double>& apdq);

void 
fflux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q);

#endif // __cold__
