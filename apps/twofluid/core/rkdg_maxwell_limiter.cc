#include "miniwarpx.h"
#include "utils.h"
#include "rkdg_algo.h"
#include "twofluid.h"
#include <math.h>

/*
  Applies characteristic based limiters to the conserved
  variables. This routine is specific to the Maxwell equations.

*/
void
rkdg_maxwell_limiter(Run_Data& rd, FArray<double>& q, double dt)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mwave = rd.mwave;
    
    double dx = rd.dx;
    double M = rd.mmM; // for use in min-mod function

    FArray<double> delta( Range(1,6),0.0 );

    // speed of light
    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;
    double c0 = tfv->c0;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    // attach wave to rws->wave
    FArray<double> wave = rws->wave;

    double a1,a2,a3,a4,a5,a6;
    double a1_1,a2_1,a3_1,a4_1,a5_1,a6_1;
    double a1_p,a2_p,a3_p,a4_p,a5_p,a6_p;
    double a1_m,a2_m,a3_m,a4_m,a5_m,a6_m;

    // compute waves
    for(int i=1; i<=mx; i++)
    {
        // first split the linear term into characteristic
        delta(1) = q(i,1,2);
        delta(2) = q(i,2,2);
        delta(3) = q(i,3,2);
        delta(4) = q(i,4,2);
        delta(5) = q(i,5,2);
        delta(6) = q(i,6,2);

        a1_1 = 0.5*(delta(3)/c0+delta(5));
        a2_1 = 0.5*(-delta(2)/c0+delta(6));
        a3_1 = delta(4);
        a4_1 = delta(1);
        a5_1 = 0.5*(-delta(3)/c0+delta(5));
        a6_1 = 0.5*(delta(2)/c0+delta(6));

        // next split the forward difference of averages into characteristic
        delta(1) = q(i+1,1,1) - q(i,1,1);
        delta(2) = q(i+1,2,1) - q(i,2,1);
        delta(3) = q(i+1,3,1) - q(i,3,1);
        delta(4) = q(i+1,4,1) - q(i,4,1);
        delta(5) = q(i+1,5,1) - q(i,5,1);
        delta(6) = q(i+1,6,1) - q(i,6,1);

        a1_p = 0.5*(delta(3)/c0+delta(5));
        a2_p = 0.5*(-delta(2)/c0+delta(6));
        a3_p = delta(4);
        a4_p = delta(1);
        a5_p = 0.5*(-delta(3)/c0+delta(5));
        a6_p = 0.5*(delta(2)/c0+delta(6));

        // next split the coefficient of the backward difference of averages
        delta(1) = q(i,1,1) - q(i-1,1,1);
        delta(2) = q(i,2,1) - q(i-1,2,1);
        delta(3) = q(i,3,1) - q(i-1,3,1);
        delta(4) = q(i,4,1) - q(i-1,4,1);
        delta(5) = q(i,5,1) - q(i-1,5,1);
        delta(6) = q(i,6,1) - q(i-1,6,1);

        a1_m = 0.5*(delta(3)/c0+delta(5));
        a2_m = 0.5*(-delta(2)/c0+delta(6));
        a3_m = delta(4);
        a4_m = delta(1);
        a5_m = 0.5*(-delta(3)/c0+delta(5));
        a6_m = 0.5*(delta(2)/c0+delta(6));

        // use minmod function to compute actual coefficients
        a1 = mminmod(a1_1,a1_p,a1_m,dx,M);
        a2 = mminmod(a2_1,a2_p,a2_m,dx,M);
        a3 = mminmod(a3_1,a3_p,a3_m,dx,M);
        a4 = mminmod(a4_1,a4_p,a4_m,dx,M);
        a5 = mminmod(a5_1,a5_p,a5_m,dx,M);
        a6 = mminmod(a6_1,a6_p,a6_m,dx,M);

        // compute waves

        // wave 1 with speed -c
        wave(i,1,1) = 0;
        wave(i,2,1) = -a2*c0;
        wave(i,3,1) = a1*c0;
        wave(i,4,1) = 0;
        wave(i,5,1) = a1;
        wave(i,6,1) = a2;

        // wave 2 with speed 0
        wave(i,1,2) = a4;
        wave(i,2,2) = 0;
        wave(i,3,2) = 0;
        wave(i,4,2) = a3;
        wave(i,5,2) = 0;
        wave(i,6,2) = 0;

        // wave 3 with speed +c
        wave(i,1,3) = 0;
        wave(i,2,3) = a6*c0;
        wave(i,3,3) = -a5*c0;
        wave(i,4,3) = 0;
        wave(i,5,3) = a5;
        wave(i,6,3) = a6;

        // use the computed waves to correct the coefficient of the
        // linear term
        for(int m=1; m<=meqn; m++)
        {
            q(i,m,2) = 0.0; // set coefficient to 0
            for(int mw=1; mw<=mwave; mw++)
                q(i,m,2) += wave(i,m,mw);
        }
    }
}
