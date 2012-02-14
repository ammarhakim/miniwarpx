#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, p,mez,ex,by,rhoe,rhoi;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = tfv->mi;
    double me = tfv->me;
    double qi = tfv->qi;
    double qe = tfv->qe;

    double rad = 0.25;
    double ratio_me_mi = me/mi;
    double alpha = 0.1;  
    double j0 = 0.5;
    double p0 = (1.0+alpha)*0.25*j0*j0*rad*rad;
    double pa = alpha*.25*j0*j0*rad*rad;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        ex = 0.0;
        
        // Inside the pinch
        if (xcell<=rad)
        {
            by = 0.5*xcell*j0;
            mez = j0*me/qe;
            p = p0-0.25*j0*j0*xcell*xcell;
            rhoe = me*p/p0;
            rhoi = mi*p/p0;
        }
        // Outside the pinch
        else
        {
            by = 0.5*rad*rad/xcell*j0;
            mez = 0.0;
            p = pa;
            rhoe = me*p/p0;
            rhoi = mi*p/p0;
        }

        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = rhoe;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = mez;
        q(i,5)  = p/gas_gamma1 + 0.5*(pow(q(i,2),2)+(q(i,3),2)
                                      +(q(i,4),2))/q(i,1);
        q(i,6)  = rhoi;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        q(i,10) = p/gas_gamma1 +  0.5*(pow(q(i,7),2)+(q(i,8),2)
                                       +(q(i,9),2))/q(i,1);
        q(i,11) = ex;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = by;
        q(i,16) = 0.0;
    }
}
