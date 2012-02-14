#include "miniwarpx.h"
#include "euler.h"

// flux function for Nozzle flow
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    double rho,u,v,w,E,p;
    double xcell, area;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    // gas constant
    double gas_gamma = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // xcell and area needed because conserved variables and 
        // fluxes have an area multiplied in them
        xcell = rd.xcoords(i);
        area = 1.398+0.347*tanh(0.8*xcell-4.); 

        // compute primitive variables
        rho = q(i,1)/area;
        u   = q(i,2)/(rho*area);
        v   = q(i,3)/(rho*area);
        w   = q(i,4)/(rho*area);
        E   = q(i,5)/area;

        p = gas_gamma1*(E-0.5*rho*(u*u+v*v+w*w));
        
        // compute flux
        fx(i,1) = rho*u*area;
        fx(i,2) = (rho*u*u + p)*area;
        fx(i,3) = rho*u*v*area;
        fx(i,4) = rho*u*w*area;
        fx(i,5) = (E+p)*u*area;
    }
}
