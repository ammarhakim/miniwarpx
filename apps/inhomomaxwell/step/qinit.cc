#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "maxwell.h"

#define is_even(x) (x%2) ? 0 : 1

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double Ey, Hz;
    IMaxwell_Vars *v = (IMaxwell_Vars*) rd.mvar;
    double thick = 0.5;

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        if (xcell < -2.5)
            Ey = 1.0;
        else
            Ey = -1.0;
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

        if(xcell < 0.0)
        {
            v->ep(i) = 1.0;
            v->mu(i) = 1.0;
        }
        else
        {
            v->ep(i) = 1.0;
            v->mu(i) = 1.0;
        }
    }
}
