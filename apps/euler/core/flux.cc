#include "miniwarpx.h"
#include "euler.h"

// flux function for Burger's equation
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    double rho,u,v,w,E,p;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    // gas constant
    double gas_gamma = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        rho = q(i,1);
        u   = q(i,2)/rho;
        v   = q(i,3)/rho;
        w   = q(i,4)/rho;
        E   = q(i,5);

        p = gas_gamma1*(E-0.5*rho*(u*u+v*v+w*w));
        
        // compute flux
        fx(i,1) = rho*u;
        fx(i,2) = rho*u*u + p;
        fx(i,3) = rho*u*v;
        fx(i,4) = rho*u*w;
        fx(i,5) = (E+p)*u;
    }
}
