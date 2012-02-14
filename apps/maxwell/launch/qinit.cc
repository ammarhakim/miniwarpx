#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "maxwell.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;

    const double pi = 3.141592654;

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);
        // electric field
        q(i,1) = 0.0;
        q(i,2) = 0.0;
        q(i,3) = 0.0;
        // magnetic field
        q(i,4) = 0.0;
        q(i,5) = 0.0;
        q(i,6) = 0.0;
    }

}
