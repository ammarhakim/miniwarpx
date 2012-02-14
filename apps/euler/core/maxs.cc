#include "miniwarpx.h"
#include "euler.h"
#include <iostream>

inline
static
double
pressure(double gas_gamma, const FArray<double>& q, int i)
{
    // rho*(u^2+v^2+w^2)
    double rhou2 = (pow(q(i,2),2) + pow(q(i,3),2) + pow(q(i,4),2))/q(i,1);

    // compute pressure
    return (gas_gamma-1)*(q(i,5) - 0.5*rhou2);
}

/*
  Computes the maxmimum wave speed (Fast Magnetosonic speed) for MHD equations

  Parameters
  ----------

  rd [in]  - Input data for simulation
  s  [out] - Maximum wave speed. s(i) is the wave speed in cell i
  q  [in]  - Conserved variable at which speed is to be computed

*/

void 
maxs(const Run_Data& rd, FArray<double>& s, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma; // gas constant

    double rho,u,v,w,cs,pr;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        rho = q(i,1);
        u = q(i,2)/rho;
        v = q(i,3)/rho;
        w = q(i,4)/rho;

        pr = pressure(gas_gamma,q,i); // fluid pressure
        cs = sqrt(gas_gamma*pr/rho); // sound speed^2

        // compute maximum wave speed, ignoring sign. Note that if u<0
        // then fabs(u-cs) > fabs(u+cs) and hence the following test
        // is needed
        if (fabs(u+cs) > fabs(u-cs))
            s(i) = u+cs;
        else
            s(i) = u-cs;
    }
}
