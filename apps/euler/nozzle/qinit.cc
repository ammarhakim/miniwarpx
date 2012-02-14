#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, area, M,c,u,p, rho,rhou,e;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma;

    // set initial conditions for density, momentum, and energy
    // using specified Mach number and pressure
    M = 1.26;
    p    = 1.0;
    rho  = 1.0;
    // speed of sound at specified conditions
    c    = sqrt(gas_gamma*p/rho);
    u    = 0.0;
    rhou = rho*u;
    e    = p/(gas_gamma-1.) + 0.5*rhou*u;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell  = xloc(i);
        area   = 1.398+0.347*tanh(0.8*xcell-4.);

        // compute the conserved variables
        q(i,1) = rho*area;
        q(i,2) = rhou*area;
        q(i,3) = 0.0;
        q(i,4) = 0.0;
        q(i,5) = e*area;
    }
}
