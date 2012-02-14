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
    double D = 0.2;

    IMaxwell_Vars *v = (IMaxwell_Vars*) rd.mvar;

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        Ey = exp(-pow(xcell-1.0,2)/(D*D));//*cos(2*pi*xcell);
        //Ey = sin(2.*pi*xcell/10.0);
        // right going wave
        Hz =  Ey;

        q(i,1) = Ey;
        q(i,2) = Hz;
    }

    double thick = 0.5;
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
            v->ep(i) = 2.0;
            v->mu(i) = 1.0;
        }
    }
    write_coeffs(rd, "epsmu");
}
