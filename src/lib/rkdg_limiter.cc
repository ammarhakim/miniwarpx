#include "miniwarpx.h"
#include "rkdg_algo.h"
#include "utils.h"
#include <math.h>

/* 
   Characteristic based limiters may be provided by the user optimized
   for the particular equation system being solved
 */
extern
void rkdg_limiter_characteristics(Run_Data& rd, FArray<double>& q, double dt);

/*
  Applies component based limiters to the conserved variables. This is
  much faster than doing charateristic based limiters, however is not
  TVDM. These limiters are useful when simple Lax fluxes are used to
  solve the RP and complete eigen-decomposition is not available.

 */
static
void
rkdg_limiter_components(Run_Data& rd, FArray<double>& q, double dt)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    double mmM = rd.mmM;
    bool changed;

    // loop over all cells applying component-wise limiters to each
    // equation
    double q1, q1p,q1m;
    double mm;
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
        {
            q1 = q(i,m,2); // coefficient of linear term
            q1p = q(i+1,m,1)-q(i,m,1); // forward difference of average
            q1m = q(i,m,1)-q(i-1,m,1); // backward difference of average

            // compute limited value of linear term
            mm = mminmod(q1,q1p,q1m,rd.dx,mmM);

            // check if we need to apply limiters
            if(mm == q(i,m,2))
                changed = false; // no, we are not applying limiter
            else
                changed = true; // yes, we are applying limiter

            q(i,m,2) = mm; // limit coefficient of linear term

            // if we applied limiters to the linear term, set all
            // higher order coefficients to 0.0
            if(changed)
                for(int cc=3; cc<=rd.sp_order; cc++)
                    q(i,m,cc) = 0.0;
        }
}

/*
  Limiter for the Runge-Kutta Discontinuous akgorithm.
 */
void 
rkdg_limiter(Run_Data& rd, FArray<double>& q, double dt)
{
    if(rd.dg_limiters == 0)
        // no limiters: do nothing
        ;
    else if(rd.dg_limiters == 1)
        // characteritics based limiters
        rkdg_limiter_characteristics(rd,q,dt);
    else if(rd.dg_limiters == 2)
        // component based limiters
        rkdg_limiter_components(rd,q,dt);
    else
    {/* bad bad very bad */ }

}
