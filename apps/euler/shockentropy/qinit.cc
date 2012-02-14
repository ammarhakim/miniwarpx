#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double rhol,ul,rhoul,pl,el;
    double ur,pr;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    double xcell,sloc = -0.8;

    // data in left state:
    rhol = 3.857143;
    ul = 2.629369;
    pl = 10.333333;
    rhoul = rhol*ul;
    el = pl/gas_gamma1 + 0.5*rhol*ul*ul;

    // data in right state:
    ur = 0.0;
    pr = 1.0;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        if(xcell < sloc)
        {
            q(i,1) = rhol;
            q(i,2) = rhoul;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            q(i,5) = el;
        }
        else
        {
            q(i,1) = 1.0 + 0.2*sin(5*PI*xcell);
            q(i,2) = q(i,1)*ur;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            q(i,5) = pr/gas_gamma1 + 0.5*q(i,1)*ur*ur;
        }
    }
}
