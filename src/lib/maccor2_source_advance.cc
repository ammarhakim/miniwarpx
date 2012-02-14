#include "miniwarpx.h"
#include "maccor_algo.h"

#include <math.h>

/*
  Solves the ODE dq/dt = s from time 'tcurr' to 'tcurr+dt'. This
  particular routine uses a 4th order Runge-Kutta method to solve the
  ODE. The user can provide his own wave_source_advance function if he
  wants to use some other method.

  Parameters
  ----------

  rd [in]     - Simulation parameters
  q  [in/out] - On input contains solution at time 'tcurr'. On output contains 
                solution at time 'tcurr+dt'.
  tcurr [in]  - Current time at which source term is being calculated
  dt    [in]  - Time step to advace ODE by

 */
void 
maccor2_source_advance(Run_Data& rd, FArray<double>& q,double tcurr, double dt)
{
    int mx = rd.mx;
    int meqn = rd.meqn;

    MACCOR2_Workspace *wws = (MACCOR2_Workspace*) rd.work; // get a hold of workspace

    // attach src, srct, and srcm to various array already allocated
    // in MACCOR2_Workspace
    FArray<double> src = wws->f;
    FArray<double> srct = wws->fl;
    FArray<double> srcm = wws->fr;

    double tstart;
    double hh,h6;

    tstart = rd.tcurrent; // need to store this to restore it while
                          // leaving this function

    hh = dt/2.0;
    h6 = dt/6.0;

    //
    // Runge-Kutta stage 1
    //
    
    // copy conserved variables into wws->q
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            wws->q(i,m) = q(i,m,1);

    rd.src(rd,src,wws->q);
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            wws->q(i,m) = q(i,m,1) + hh*src(i,m);

    //
    // Runge-Kutta stage 2
    //
    
    rd.tcurrent = tstart + hh;
    rd.src(rd,srct,wws->q);
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            wws->q(i,m) = q(i,m,1) + hh*srct(i,m);

    //
    // Runge-Kutta stage 3
    // 

    rd.src(rd,srcm,wws->q);
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
        {
            wws->q(i,m) = q(i,m,1) + dt*srcm(i,m);
            srcm(i,m) = srct(i,m) + srcm(i,m);
        }

    //
    // Runge-Kutta stage 4
    // 
    rd.tcurrent = tstart + dt;
    rd.src(rd,srct,wws->q);
    // perform final update
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            q(i,m,1) = q(i,m,1) + h6*(src(i,m)+srct(i,m)+2.*srcm(i,m));

    rd.tcurrent = tstart; // restore rd.tcurrent

}
