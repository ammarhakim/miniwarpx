#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "cold.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, ne, ue, ex;
    const double D_PI = 3.141592654;
    double n0 = 4.0;
    double w, k;

    w = sqrt(n0);
    k = 4*D_PI;

    for(int i=1; i<=rd.mx; ++i)
    {
        xcell = xloc(i);
        ue = 1e-4*sin(k*xcell);
        ne = -n0*ue*(k/w);
        ex = 1e-4*n0/w*cos(k*xcell);

        // set electron variables        
        q(i,1) = n0 + ne;
        q(i,2) = q(i,1)*ue;
        q(i,3) = 0.0;
        q(i,4) = 0.0;
        // set EM variables
        q(i,5) = ex;
        q(i,6) = 0.0;
        q(i,7) = 0.0;
        q(i,8) = 0.0;
        q(i,9) = 0.0;
        q(i,10) = 0.0;
    }
}
