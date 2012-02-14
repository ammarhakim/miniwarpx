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

    double n0 = 1.0;
    double p0 = 1.0;
    double alpha = 0.1;   // this is 1/(rad^2)
    double beta = 0.1;    // this is pback/p0

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);

        p  = p0/exp(pow(xcell,2)*alpha) + p0*beta;

        ni = n0/exp(pow(xcell,2)*alpha) + n0*beta;

        ne = (pow(n0,2)*pow(qi,2)*pow(1 + exp(pow(xcell,2)*alpha)*beta,3)
              - 2*exp(pow(xcell,2)*alpha)*p0*alpha*(-1 + exp(pow(xcell,2)
              *alpha)*(-1 + pow(xcell,2)*alpha)*beta)*epsilon0)
            /(exp(pow(xcell,2)*alpha)*n0*pow(
                  qi + exp(pow(xcell,2)*alpha)*qi*beta,2));


         mez= -((2*p0*xcell*alpha*(pow(n0,2)*pow(qi,2)*pow(1 
               + exp(pow(xcell,2)*alpha)*beta,3) + 
               exp(pow(xcell,2)*alpha)*p0*alpha*(1 + exp(pow(xcell,2)*alpha)
               *(beta - pow(xcell,2)*alpha*beta))*epsilon0))/(exp(pow(xcell,2)*alpha)
               *pow(n0,2)*pow(qi + exp(pow(xcell,2)*alpha)*qi*beta,3)*
               sqrt((p0*(2 - (2*(1 + pow(xcell,2)*alpha))/exp(pow(xcell,2)*alpha) + 
               (p0*pow(xcell,4)*pow(alpha,3)*epsilon0)/(pow(n0,2)*pow(qi,2)
               *pow(1 + exp(pow(xcell,2)*alpha)*beta,2))))/(pow(xcell,2)*alpha*1))
                                                                *pow(1,2)));

        by = sqrt((2.*p0)/(pow(xcell,2)*alpha*1.) + (p0*alpha*(-2.*pow(n0,2)*pow(qi,2)*
              (1. + pow(xcell,2)*alpha)*pow(1. + exp(pow(xcell,2)*alpha)*beta,2) + 
              exp(pow(xcell,2)*alpha)*p0*pow(xcell,4)*pow(alpha,3)*epsilon0))/
              (exp(pow(xcell,2)*alpha)*pow(n0,2)*pow(qi,2)*pow(xcell,2)*
              pow(alpha + exp(pow(xcell,2)*alpha)*alpha*beta,2)*1.));

        ex = -((p0*xcell*alpha)/(exp(pow(xcell,2)*alpha)*qi
                                 *(n0/exp(pow(xcell,2)*alpha) + n0*beta)));

        pe = 0.5*p;

        pi = pe;

        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = me*ne;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = mez;
        q(i,5)  = pe/gas_gamma1 + 0.5*(pow(q(i,2),2)+pow(q(i,3),2)
                                       +pow(q(i,4),2))/q(i,1);
        q(i,6)  = mi*ni;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        q(i,10) = pi/gas_gamma1 +  0.5*(pow(q(i,7),2)+pow(q(i,8),2)
                                       +pow(q(i,9),2))/q(i,6);
        q(i,11) = ex;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = by;
        q(i,16) = 0.0;
    }
}
