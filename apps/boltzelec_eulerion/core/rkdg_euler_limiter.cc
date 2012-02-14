#include "miniwarpx.h"
#include "utils.h"
#include "rkdg_algo.h"
#include "euler.h"
#include <math.h>

/*
  Applies characteristic based limiters to the conserved
  variables. This routine is specific to the Euler equations.

 */
void
rkdg_limiter_characteristics(Run_Data& rd, FArray<double>& q, double dt)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mbc = rd.mbc;
    //int mwave = rd.mwave;
    int mwave = 3;
    
    double dx = rd.dx;
    double M = rd.mmM; // for use in min-mod function

    // attach qt to the cell averages
    FArray<double> qt(Range(1-mbc,mx+mbc, 1,meqn), q.data());
    FArray<double> wave(Range(1,meqn, 1,mwave),0.0);

    double rho,pr;
    double a1,a2,a3,a4,a5;
    double a1_1,a2_1,a3_1,a4_1,a5_1;
    double a1_p,a2_p,a3_p,a4_p,a5_p;
    double a1_m,a2_m,a3_m,a4_m,a5_m;

    FArray<double> delta(Range(1,5),0.0);

    // gas constant
    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    // gas constant
    double gas_gamma = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;

    // compute quantities useful for characteristic splitting
    for(int i=2-mbc; i<=mx+mbc; i++)
    {

        // density
        rho = qt(i,1);

        // velocity components
        ev->u(i) = qt(i,2)/rho;
        ev->v(i) = qt(i,3)/rho;
        ev->w(i) = qt(i,4)/rho;

        // pressure
        pr = gas_gamma1*(qt(i,5) 
                         - 0.5*(pow(qt(i,2),2) +
                                pow(qt(i,3),2) +
                                pow(qt(i,4),2))/rho);
        
        ev->u2v2w2(i) = pow(ev->u(i),2) + pow(ev->v(i),2) + pow(ev->w(i),2);

        // enthalpy
        ev->enth(i) = (qt(i,5)+pr)/rho;

        // speed of sound
        a2 = gas_gamma1*(ev->enth(i) - .5*ev->u2v2w2(i));
        ev->a(i) = sqrt(a2);

        ev->g1a2(i) = gas_gamma1 / a2;
        ev->euv(i)  = ev->enth(i) - ev->u2v2w2(i);
    }

    // compute waves
    for(int i=1; i<=mx; i++)
    {
        // compute coefficients of the 5 eigenvectors

        // first split the coefficient of the linear term
        delta(1) = q(i,1,2);
        delta(2) = q(i,2,2);
        delta(3) = q(i,3,2);
        delta(4) = q(i,4,2);
        delta(5) = q(i,5,2);

        a4_1 = ev->g1a2(i) * (ev->euv(i)*delta(1)
                              + ev->u(i)*delta(2) 
                              + ev->v(i)*delta(3) 
                              + ev->w(i)*delta(4)
                              - delta(5));
        a2_1 = delta(3) - ev->v(i)*delta(1);
        a3_1 = delta(4) - ev->w(i)*delta(1);
        a5_1 = (delta(2) + (ev->a(i)-ev->u(i))*delta(1) - ev->a(i)*a4_1) / (2.0*ev->a(i));
        a1_1 = delta(1) - a4_1 - a5_1;

        // next split the coefficient of the forward difference of averages
        delta(1) = q(i+1,1,1) - q(i,1,1);
        delta(2) = q(i+1,2,1) - q(i,2,1);
        delta(3) = q(i+1,3,1) - q(i,3,1);
        delta(4) = q(i+1,4,1) - q(i,4,1);
        delta(5) = q(i+1,5,1) - q(i,5,1);

        a4_p = ev->g1a2(i) * (ev->euv(i)*delta(1)
                              + ev->u(i)*delta(2) 
                              + ev->v(i)*delta(3) 
                              + ev->w(i)*delta(4)
                              - delta(5));
        a2_p = delta(3) - ev->v(i)*delta(1);
        a3_p = delta(4) - ev->w(i)*delta(1);
        a5_p = (delta(2) + (ev->a(i)-ev->u(i))*delta(1) - ev->a(i)*a4_p) / (2.0*ev->a(i));
        a1_p = delta(1) - a4_p - a5_p;

        // next split the coefficient of the backward difference of averages
        delta(1) = q(i,1,1) - q(i-1,1,1);
        delta(2) = q(i,2,1) - q(i-1,2,1);
        delta(3) = q(i,3,1) - q(i-1,3,1);
        delta(4) = q(i,4,1) - q(i-1,4,1);
        delta(5) = q(i,5,1) - q(i-1,5,1);

        a4_m = ev->g1a2(i) * (ev->euv(i)*delta(1)
                              + ev->u(i)*delta(2) 
                              + ev->v(i)*delta(3) 
                              + ev->w(i)*delta(4)
                              - delta(5));
        a2_m = delta(3) - ev->v(i)*delta(1);
        a3_m = delta(4) - ev->w(i)*delta(1);
        a5_m = (delta(2) + (ev->a(i)-ev->u(i))*delta(1) - ev->a(i)*a4_m) / (2.0*ev->a(i));
        a1_m = delta(1) - a4_m - a5_m;

        // use minmod function to compute actual coefficients
        a1 = mminmod(a1_1,a1_p,a1_m,dx,M);
        a2 = mminmod(a2_1,a2_p,a2_m,dx,M);
        a3 = mminmod(a3_1,a3_p,a3_m,dx,M);
        a4 = mminmod(a4_1,a4_p,a4_m,dx,M);
        a5 = mminmod(a5_1,a5_p,a5_m,dx,M);

        // compute the waves

        // Wave 1: eigenvalue is u-c
        wave(1,1) = a1;
        wave(2,1) = a1*(ev->u(i)-ev->a(i));
        wave(3,1) = a1*ev->v(i);
        wave(4,1) = a1*ev->w(i);
        wave(5,1) = a1*(ev->enth(i) - ev->u(i)*ev->a(i));

        // Wave 2: the 3 eigenvectors corresponding to the repeated
        // eigenvalue u,u,u are lumped together into a single wave
        wave(1,2) = a4;
        wave(2,2) = a4*ev->u(i);
        wave(3,2) = a4*ev->v(i)	 	 + a2;
        wave(4,2) = a4*ev->w(i)	 	 + a3;
        wave(5,2) = a4*0.5*ev->u2v2w2(i)  + a2*ev->v(i) + a3*ev->w(i);

        // Wave 3: eigenvalue is u+c
        wave(1,3) = a5;
        wave(2,3) = a5*(ev->u(i)+ev->a(i));
        wave(3,3) = a5*ev->v(i);
        wave(4,3) = a5*ev->w(i);
        wave(5,3) = a5*(ev->enth(i)+ev->u(i)*ev->a(i));

        // use the computed waves to correct the coefficient of the
        // linear term
        for(int m=1; m<=meqn; m++)
        {
            q(i,m,2) = 0.0; // set coefficient to 0
            for(int mw=1; mw<=mwave; mw++)
                q(i,m,2) += wave(m,mw);
        }
    }
}
