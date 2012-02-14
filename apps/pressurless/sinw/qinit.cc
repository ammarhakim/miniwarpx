#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double pi = 3.1415292654;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        q(i,1) = 1.0;
        q(i,2) = 0.1*sin(2*pi*xcell);
    }
}
