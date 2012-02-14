#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "maxwell.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double Ey,Ez,By,Bz;

    const double pi = 3.141592654;
    double D = 0.1;

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        Ey = exp(-pow(xcell-0.5,2)/(D*D))*cos(2*pi*xcell);
        Ez = exp(-pow(xcell-0.5,2)/(D*D))*sin(2*pi*xcell);

        // right going wave
        By = -Ez;
        Bz =  Ey;

        // electric field
        q(i,1) = 0.0;
        q(i,2) = Ey;
        q(i,3) = Ez;
        // magnetic field
        q(i,4) = 0.0;
        q(i,5) = By;
        q(i,6) = Bz;
    }

}
