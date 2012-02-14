#include "miniwarpx.h"
#include "utils.h"
#include "rkdg_algo.h"
#include "twofluid.h"
#include <math.h>

/*
  Applies characteristic based limiters to the conserved
  variables. This routine is specific to the Two-Fluid equations. The
  limiters are applied seprately to each fluid and EM fields by
  calling the appropriate limiter function for each equation system.

 */
void
rkdg_limiter_characteristics(Run_Data& rd, FArray<double>& q, double dt)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    int nc = rd.sp_order;

    int start[3]; // for use in indicating starting indices of sliced arraysa
    FArray<double> q_s; // slice of q to pass to limiters

    start[0] = 1-mbc;
    start[1] = 1;
    start[2] = 1;

    // slice array to pass to Euler equation limiter
    q_s = q.slice( Range(1-mbc,mx+mbc, 1,5, 1,nc), start);
    rd.meqn = 5;
    rd.mwave = 3;

    //
    // call limiter for Euler equations for electron fluid
    //
    rkdg_euler_limiter(rd,q_s,dt);

    // slice array to pass to Euler equation limiter
    q_s = q.slice( Range(1-mbc,mx+mbc, 6,10, 1,nc), start);
    rd.meqn = 5;
    rd.mwave = 3;

    //
    // call limiter for Euler equations for ion fluid
    //
    rkdg_euler_limiter(rd,q_s,dt);

    // slice array to pass to Maxwell equation limiter
    q_s = q.slice( Range(1-mbc,mx+mbc, 11,16, 1,nc), start);
    rd.meqn = 6;
    rd.mwave = 3;

    //
    // call limiter for Maxwell equations
    //
    rkdg_maxwell_limiter(rd,q_s,dt);

    // restore values of meqn and mwave
    rd.meqn = 16;
    rd.mwave = 9;
}
