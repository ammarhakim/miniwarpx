#ifndef __hallmhd_h__
#define __hallmhd_h__

#include "farray.h"

// data-structure to hold variables for Euler equation application
struct Hallmhd_Vars
{
    double gas_gamma; // adiabatic index of gas
    double c0;        // speed of light
    double qi;        // ion charge
    double qe;        // electron charge
    double mi;        // ion mass
    double me;        // electron mass
    double epsilon0; 
    int is_radial;    // for radial component
    double ur;        // reference velocity (thermal vel.)
    double betar;     // reference beta
    double rli;       // reference ion larmor radius

    // workspace needed by the hallmhd equation Reimann solvers
    FArray<double> u2v2w2,u,v,w,enth,a,g1a2,euv;
};

// void  rp_maxwell(Run_Data&, FArray<double>&, FArray<double>&, FArray<double>&,
//                  FArray<double>&, FArray<double>& ,
//                  FArray<double>&, FArray<double>&);

 void  rp_euler(Run_Data&, FArray<double>&, FArray<double>&, FArray<double>&,
                FArray<double>&, FArray<double>&,
                FArray<double>&, FArray<double>&);

// void rkdg_euler_limiter(Run_Data&, FArray<double>&, double);
// void rkdg_maxwell_limiter(Run_Data&, FArray<double>&, double);


#endif // __hallmhd_h__
