#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "maxwell.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double Ey,Hz;
    const double pi = 3.141592654;
    double D = 0.25;

    IMaxwell_Vars *v = (IMaxwell_Vars*) rd.mvar;

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        Ey = exp(-pow(xcell-0.5,2)/(D*D));//*cos(2*pi*xcell);
        // right going wave
        Hz =  Ey;

        q(i,1) = Ey;
        q(i,2) = Hz;
    }

    // set epsilon and mu
    for(int i=1-rd.mbc; i<=rd.mx+rd.mbc; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        if(xcell < 3.0)
        {
            v->ep(i) = 1.0;
            v->mu(i) = 1.0;
        }
        else if ((xcell > 3.0) && (xcell < 6.0))
        {
            v->ep(i) = 2.0;
            v->mu(i) = 1.0;
        }
        else
        {
            v->ep(i) = 1.0;
            v->mu(i) = 1.0;
        }
    }

}
