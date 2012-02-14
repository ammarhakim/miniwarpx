#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double pi2 = 2.0*3.141592654;
    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        q(i,1) = 0.8*sin(pi2*xcell);
    }
}
