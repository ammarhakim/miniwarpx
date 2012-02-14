#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double r0 = 0.5, u0;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        if(xcell < -0.5) 
            u0 = -0.5;
        else if((-0.5 < xcell) && (xcell < 0))
            u0 = 0.4;
        else if ((0 < xcell) && (xcell < 0.8))
            u0 = 0.4-xcell;
        else
            u0 = -0.4;

        q(i,1) = r0;
        q(i,2) = r0*u0;
    }
}
