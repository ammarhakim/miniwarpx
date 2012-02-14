#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "mhd.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double sloc,pr;

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma; // gas constant
    double gas_gamma1 = gas_gamma - 1;

    sloc = 0.0;

    // Brio-Wu MHD shock-tube problem
    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        if (xcell < sloc)
        {
            q(i,1) = 1.0;
            q(i,2) = 0.0;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            
            q(i,6) = 0.75;
            q(i,7) = 1.0;
            q(i,8) = 0.0;
            
            pr = 1.0;
            q(i,5) = pr/gas_gamma1 + 0.5*(pow(q(i,6),2)+pow(q(i,7),2)+pow(q(i,8),2));
        }
        else
        {
            q(i,1) = 0.125;
            q(i,2) = 0.0;
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            
            q(i,6) = 0.75;
            q(i,7) = -1.0;
            q(i,8) = 0.0;
            
            pr = 0.1;
            q(i,5) = pr/gas_gamma1 + 0.5*(pow(q(i,6),2)+pow(q(i,7),2)+pow(q(i,8),2));
        }
    }
}
