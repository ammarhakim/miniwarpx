#include "miniwarpx.h"
#include "utils.h"
#include "rkdg_algo.h"
#include "strainwave.h"
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
    int mwave = rd.mwave;
    
    double dx = rd.dx;
    double M = rd.mmM; // for use in min-mod function

    double a1_1,a2_1,a1_p,a2_p,a1_m,a2_m;
    double a1,a2,c0,Z,dsig_dstrain;
    Strainwave_Vars *sv = (Strainwave_Vars*) rd.mvar;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    // attach qt to the cell averages
    FArray<double> qt(Range(1-mbc,mx+mbc, 1,meqn), q.data());
    // attach wave to rws->wave
    FArray<double> wave = rws->wave;

    FArray<double> delta(Range(1,2),0.0);

    // compute quantities
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        // for linear problems beta=0, so can use these statements
        // for both linear and nonlinear problems
        dsig_dstrain = sv->modul(i) 
            + 2*sv->beta*sv->modul(i)*sv->modul(i)*qt(i,1);
        c0 = sqrt(dsig_dstrain/sv->rho(i));
        Z  = sqrt(sv->rho(i)*dsig_dstrain);

        // first split the coefficient of the linear term
        delta(1) = q(i,1,2);
        delta(2) = q(i,2,2);

        a1_1 = 0.5*(delta(1)+delta(2)/Z);
        a2_1 = 0.5*(delta(1)-delta(2)/Z);

        // next split the coefficient of the forward difference of averages
        delta(1) = q(i+1,1,1) - q(i,1,1);
        delta(2) = q(i+1,2,1) - q(i,2,1);

        a1_p = 0.5*(delta(1)+delta(2)/Z);
        a2_p = 0.5*(delta(1)-delta(2)/Z);

        // next split the coefficient of the backward difference of averages
        delta(1) = q(i,1,1) - q(i-1,1,1);
        delta(2) = q(i,2,1) - q(i-1,2,1);

        a1_m = 0.5*(delta(1)+delta(2)/Z);
        a2_m = 0.5*(delta(1)-delta(2)/Z);

        // use minmod function to compute actual coefficients
        a1 = mminmod(a1_1,a1_p,a1_m,dx,M);
        a2 = mminmod(a2_1,a2_p,a2_m,dx,M);

        // compute the waves
        // wave 1 with speed -c
        wave(i,1,1) = a1;
        wave(i,2,1) = a1*Z;

        // wave 2 with speed +c
        wave(i,1,2) = a2;
        wave(i,2,2) = -a2*Z;

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
