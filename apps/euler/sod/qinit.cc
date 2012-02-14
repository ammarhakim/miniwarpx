#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "euler.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double sloc = 0.5;
    double rhol, rhoul, pl, el;
    double rhor, rhour, pr, er;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma = ev->gas_gamma;

    // set right and left initial states
    rhol = 3.0;
    rhoul = 0.0;
    pl = 3.0;
    el = pl/(gas_gamma-1) + 0.;

    rhor = 1.0;
    rhour = 0.0;
    pr = 1.0;
    er = pr/(gas_gamma-1) + 0.;

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
            q(i,1) = rhor;
            q(i,2) = rhour;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            q(i,5) = er;
        }
    }
}
