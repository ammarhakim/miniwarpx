#include "miniwarpx.h"
#include "mhd.h"
#include <iostream>

inline
static
double
pressure(double gas_gamma, const FArray<double>& q, int i)
{
    // rho*(u^2+v^2+w^2)
    double rhou2 = (pow(q(i,2),2) + pow(q(i,3),2) + pow(q(i,4),2))/q(i,1);
    // Bx^2+By^2+Bz^2
    double B2 = pow(q(i,6),2) + pow(q(i,7),2) + pow(q(i,8),2);

    // compute pressure
    return (gas_gamma-1)*(q(i,5) - 0.5*rhou2 - 0.5*B2);
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

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma; // gas constant

    double rho,u,v,w,pr,bx,by,bz;
    double b2,cs2,ca,cfs;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        rho = q(i,1);
        u = q(i,2)/rho;
        v = q(i,3)/rho;
        w = q(i,4)/rho;

        pr = pressure(gas_gamma,q,i); // fluid pressure

        bx = q(i,6);
        by = q(i,7);
        bz = q(i,8);

        b2 = bx*bx + by*by + bz*bz;

        // compute maximum wave speed (fast magnetosonic wave)

        cs2 = gas_gamma*pr/rho; // sound speed^2
        ca = bx*bx/rho; // Alfven speed^2

        // fast magnetosonic speed
        cfs = sqrt(0.5*(b2/rho + cs2 + sqrt(pow(b2/rho+cs2, 2.0) - 4.0*cs2*ca)));

        // compute maximum wave speed, ignoring sign. Note that if u<0
        // then fabs(u-cfs) > fabs(u+cfs) and hence the following test
        // is needed
        if (fabs(u+cfs) > fabs(u-cfs))
            s(i) = u+cfs;
        else
            s(i) = u-cfs;
    }
}
