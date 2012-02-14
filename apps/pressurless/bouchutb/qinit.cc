#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        if(xcell < 0.0) 
        {
            q(i,1) = 1e-20; // rho = 0.0
            q(i,2) = 0.0;
        }
        else
        {
            q(i,1) = 1.0;
            q(i,2) = 0.0;
        }
    }
}
