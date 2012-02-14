#ifndef __euler_h__
#define __euler_h__

#include "farray.h"

// data-structure to hold variables for Euler equation application
struct Euler_Vars
{
    double gas_gamma; // adiabatic index of gas
    double qbym; // charge to mass ratio

    int is_radial; // = 1 include radial source terms

    // workspace needed by the euler equation Reimann solvers
    FArray<double> u2v2w2,u,v,w,enth,a,g1a2,euv;

    // worspace to store "magnetic" field
    FArray<double> bf; 
};

struct GSL_Data 
{
    double *bfield;
    double qbym;
};

#endif // __euler_h__
