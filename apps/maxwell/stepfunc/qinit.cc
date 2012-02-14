#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "maxwell.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double sloc = 0.5;
    double Ex, Ey, Ez, Bx, By, Bzl, Bzr;

    Maxwell_Vars *mv = (Maxwell_Vars*) rd.mvar;

    Ex = 0;
    Ey = 0;
    Ez = 0;
    Bx = 0;
    By = 0;
    Bzl = -1;
    Bzr = 1;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        q(i,1) = Ex;
        q(i,2) = Ey;
        q(i,3) = Ez;
        q(i,4) = Bx;
        q(i,5) = By;

        if(xcell < sloc)
        {
            q(i,6) = Bzl;
        }
        else
        {
            q(i,6) = Bzr;
        }
    }
}
