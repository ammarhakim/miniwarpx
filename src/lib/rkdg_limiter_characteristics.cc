#include "miniwarpx.h"
#include "rkdg_algo.h"
#include "utils.h"
#include <math.h>


/*
  Applies characteristic based limiters to the conserved
  variables. This limiter calls the RP solver 3 times and hence and
  slow down the code considerably. It is recommended that the limiter
  be re-written for each equation system by hand to make it more
  efficient.

 */
void
rkdg_limiter_characteristics(Run_Data& rd, FArray<double>& q, double dt)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mbc = rd.mbc;
    int mwave = rd.mwave;

    double dx = rd.dx;
    double mm;
    double mmM = rd.mmM;

    bool changed = false;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    // attach qt to the cell averages
    FArray<double> qt(Range(1-mbc,mx+mbc, 1,meqn), q.data());

    // copy q1 to rws->f
    for(int i=1-mbc; i<=mx+mbc; i++)
        for(int m=1; m<=meqn; m++)
            rws->f(i,m) = q(i,m,2);

    // compute forward and backward differences of cell averages for
    // each cell
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
        {
            rws->qr(i,m) = qt(i+1,m)-qt(i,m);
            rws->ql(i,m) = qt(i,m)-qt(i-1,m);
        }

    // split q1 into waves
    rd.rp(rd,qt,qt,rws->f,rws->wave,rws->s,rws->amdq,rws->apdq);
    // split forward differences into waves
    rd.rp(rd,qt,qt,rws->ql,rws->wave1,rws->s,rws->amdq,rws->apdq);
    // split backward differences into waves
    rd.rp(rd,qt,qt,rws->qr,rws->wave2,rws->s,rws->amdq,rws->apdq);

    // now apply modified min-mod limiter to each component of each
    // wave and reconstruct cell averages
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
        {
            q(i,m,2) = 0.0; // set coefficient to 0

            changed = false;
            for(int mw=1; mw<=mwave; mw++)
            {
                // limit wave components, accumulating it to cell
                // average
                mm = mminmod(rws->wave (i,m,mw),
                             rws->wave1(i,m,mw),
                             rws->wave2(i,m,mw),
                             dx,
                             mmM);
                q(i,m,2) += mm;
                
                if(mm == rws->wave(i,m,mw))
                    changed = false || changed;
                else
                    changed = true;
            }
            
            // if we applied limiters to the linear term, set all
            // higher order coefficients to 0.0
            if(changed)
                for(int cc=3; cc<=rd.sp_order; cc++)
                    q(i,m,cc) = 0.0;
        }
}
