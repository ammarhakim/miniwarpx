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
        if((-2<xcell) && (xcell<-1))
        {
            q(i,1) = 2.0; // rho = 2
            q(i,2) = q(i,1)*1.0; // u = 1
        }
        else if((1<xcell) && (xcell<5))
        {
            q(i,1) = 1.0; // rho = 1
            q(i,2) = -q(i,1)*1.0; // u = -1
        }
        else
        {
            q(i,1) = 1e-20; // rho = 0 (note small value)
            q(i,2) = 0.0; // u = 0
        }
    }
}
