#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "cold.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, sloc;
    sloc = 2.5;

    for(int i=1; i<=rd.mx; ++i)
    {
        xcell = xloc(i);
        if(xcell < sloc)
        {
            // electron variables
            q(i,1) = 0.1;
            q(i,2) = 0.0;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            // set EM variables
            q(i,5) = 0.0;
            q(i,6) = 0.0;
            q(i,7) = 0.0;
            q(i,8) = 0.75e-2;
            q(i,9) = 0.0;
            q(i,10) = -1.0e-2;
        }
        else
        {
            // electron variables
            q(i,1) = 0.125*0.1;
            q(i,2) = 0.0;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            // set EM variables
            q(i,5) = 0.0;
            q(i,6) = 0.0;
            q(i,7) = 0.0;
            q(i,8) = 0.75e-2;
            q(i,9) = 0.0;
            q(i,10) = 1.0e-2;
        }
    }
}
