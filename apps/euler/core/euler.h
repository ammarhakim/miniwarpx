#ifndef __euler_h__
#define __euler_h__

#include "farray.h"

// data-structure to hold variables for Euler equation application
struct Euler_Vars
{
    double gas_gamma; // adiabatic index of gas

    // workspace needed by the euler equation Reimann solvers
    FArray<double> u2v2w2,u,v,w,enth,a,g1a2,euv;
};

#endif // __euler_h__
