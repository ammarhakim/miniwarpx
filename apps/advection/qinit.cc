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
        q(i,1) = exp(-10.0*pow(xcell+1,2));
    }
}
