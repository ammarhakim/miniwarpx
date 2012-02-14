#include "miniwarpx.h"
#include "rkdg_algo.h"
#include <math.h>

/**
   Computes the conserved variable at quadrature point 'c' in each
   cell.
 */
void
rkdg_eval_expansion(const Run_Data& rd, FArray<double>& qe, const FArray<double>& qi, int c)
{
    int mbc = rd.mbc;
    int mx = rd.mx;
    int meqn = rd.meqn;
    int order = rd.sp_order;

    // get a hold of the RKDG workspace 
    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    // loop over each grid point, computing 'qe' at local coordinate 'eta'
    for(int i=1-mbc; i<=mx+mbc; i++)
    {
        for(int m=1; m<=meqn; m++)
        {
            qe(i,m) = qi(i,m,1);
            for(int cc=2; cc<=order; cc++)
                qe(i,m) += qi(i,m,cc)*rws->pmx(cc,c);
        }
    }

}


/**
   Computes the conserved variable at right edge of each cell
   cell.
 */
void
rkdg_eval_expansion_right_edge(const Run_Data& rd, FArray<double>& qe, const FArray<double>& qi)
{
    int mbc = rd.mbc;
    int mx = rd.mx;
    int meqn = rd.meqn;
    int order = rd.sp_order;

    // loop over each grid point, computing 'qe' at local coordinate 'eta'
    for(int i=1-mbc; i<=mx+mbc; i++)
    {
        for(int m=1; m<=meqn; m++)
        {
            qe(i,m) = qi(i,m,1);
            for(int cc=2; cc<=order; cc++)
                qe(i,m) += qi(i,m,cc);
        }
    }
}


/**
   Computes the conserved variable at left edge of each cell
   cell.
 */
void
rkdg_eval_expansion_left_edge(const Run_Data& rd, FArray<double>& qe, const FArray<double>& qi)
{
    int mbc = rd.mbc;
    int mx = rd.mx;
    int meqn = rd.meqn;
    int order = rd.sp_order;
    int sgn;

    // loop over each grid point, computing 'qe' at local coordinate 'eta'
    for(int i=1-mbc; i<=mx+mbc; i++)
    {
        for(int m=1; m<=meqn; m++)
        {
            qe(i,m) = qi(i,m,1);
            sgn = -1;
            for(int cc=2; cc<=order; cc++)
            {
                qe(i,m) += sgn*qi(i,m,cc);
                sgn = -1*sgn;
            }
        }
    }
}
