#include <math.h>
#include "miniwarpx.h"
#include "mhd.h"

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

// flux function for MHD equations
void
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma; // gas constant

    double rho,u,v,w,pr,bx,by,bz;

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

        pr = pr + 0.5*(bx*bx + by*by + bz*bz); // total pressure

        // compute flux
        fx(i,1) = rho*u;
        fx(i,2) = rho*u*u - bx*bx + pr;
        fx(i,3) = rho*u*v - bx*by;
        fx(i,4) = rho*u*w - bx*bz;
        fx(i,5) = u*(q(i,5)+pr) - bx*(u*bx + v*by + w*bz);
        fx(i,6) = 0.0;
        fx(i,7) = u*by - v*bx;
        fx(i,8) = u*bz - w*bx;
    }
}
