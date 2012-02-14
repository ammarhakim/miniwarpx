#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell, p,mez,ne,ni,pe,pi,ex,by;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double epsilon0 = tfv->epsilon0;
    double mi = tfv->mi;
    double me = tfv->me;
    double qi = tfv->qi;

    double ratio_me_mi = me/mi;
    double n0 = 10.0;
    double p0 = 1.0;
    double alpha = 1.0;
    double beta = 0.1;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        p = p0/exp(xcell*xcell*alpha)+p0*beta;

        ne = (n0*n0*qi*qi*pow(1+exp(xcell*xcell*alpha)*beta,3)
              -2*exp(xcell*xcell*alpha)*p0*alpha
              *(-1+exp(xcell*xcell*alpha)
              *(-1+xcell*xcell*alpha)*beta)*epsilon0)
            /(exp(xcell*xcell*alpha)*n0*pow(qi+exp(xcell*xcell*alpha)*qi*beta,2));
        mez= (2*p0*xcell*alpha*(n0*n0*qi*qi*pow(1+exp(xcell*xcell*alpha)*beta,3)
              +exp(xcell*xcell*alpha)*p0*alpha*(1+exp(xcell*xcell*alpha)
              *(beta-xcell*xcell*alpha*beta))*epsilon0))
            /(exp(xcell*xcell*alpha)*n0*n0*pow(qi+exp(xcell*xcell*alpha)*qi*beta,3)
              *sqrt((p0*(2-2*(1+xcell*xcell*alpha)/exp(xcell*xcell*alpha)
              +p0*pow(xcell,4)*pow(alpha,3)*epsilon0
              /(n0*n0*qi*qi*pow(1+exp(xcell*xcell*alpha)*beta,2))))/(xcell*xcell*alpha)));
        pe = 0.5*p;

        ni = n0/exp(xcell*xcell*alpha)+n0*beta;
        pi = pe;

        ex = -p0*xcell*alpha/(exp(xcell*xcell*alpha)*qi*(n0/exp(xcell*xcell*alpha)+n0*beta));
        by = -sqrt(2*p0/(xcell*xcell*alpha)+(p0*alpha*(-2*n0*n0*qi*qi*(1+xcell*xcell*alpha)
               *pow(1+exp(xcell*xcell*alpha)*beta,2)+exp(xcell*xcell*alpha)
               *p0*pow(xcell,4)*pow(alpha,3)*epsilon0))
            /(exp(xcell*xcell*alpha)*n0*n0*qi*qi*xcell*xcell
              *pow(alpha+exp(xcell*xcell*alpha)*alpha*beta,2)));

        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = ne;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = mez;
        q(i,5)  = pe/gas_gamma1;
        q(i,6)  = ni;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        q(i,10) = pi/gas_gamma1;
        q(i,11) = ex;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = by;
        q(i,16) = 0.0;
    }
}
