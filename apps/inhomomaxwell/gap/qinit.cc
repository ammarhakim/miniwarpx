#include <stdio.h>
#include <math.h>
#include <iostream>
#include "miniwarpx.h"
#include "maxwell.h"

#define is_even(x) (x%2) ? 0 : 1

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double Ey,Hz;
    const double pi = 3.141592654;
    double thick = 0.5;

    IMaxwell_Vars *v = (IMaxwell_Vars*) rd.mvar;
    double lambda = 4.0*thick;

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        if(xcell < 5.5)
            Ey = sin(2*pi*xcell/lambda);
        else
            Ey = 0.0;
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

        if ( is_even(int(xcell/thick)) )
        {
            v->ep(i) = 1.0;
            v->mu(i) = 1.0;
        }
        else
        {
            v->ep(i) = 4.0;
            v->mu(i) = 1.0;
        }
        if (xcell < 5.5)
        {
            v->ep(i) = 1.0;
            v->mu(i) = 1.0;
        }
    }
    write_coeffs(rd, "epsmu");
}
