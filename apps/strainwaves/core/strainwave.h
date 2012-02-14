#ifndef __strainwave_h__
#define __strainwave_h__

#include "farray.h"

// data-structure to hold variables for Euler equation application
struct Strainwave_Vars
{
    double ubar, beta;
    FArray<double> modul,rho;
};

double stress(double bet,double modulus,double strain);

#endif // __strainwave_h__
