#include "miniwarpx.h"
#include "maxwell.h"

//flux function for Maxwell's equations
void
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    double Bx,By,Bz,Ex,Ey,Ez;

    Maxwell_Vars *mv = (Maxwell_Vars*) rd.mvar;
    // speed of light
    double c0 = mv->c0;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        Ex = q(i,1);
        Ey = q(i,2);
        Ez = q(i,3);
        Bx = q(i,4);
        By = q(i,5);
        Bz = q(i,6);
        
        // compute flux
        fx(i,1) = 0;
        fx(i,2) = c0*c0*Bz;
        fx(i,3) = -c0*c0*By;
        fx(i,4) = 0;
        fx(i,5) = -Ez;
        fx(i,6) = Ey;

    }
}
