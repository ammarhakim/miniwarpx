#include "miniwarpx.h"
#include "rkdg_algo.h"

/**
   Set initial conditions for the RKDG algorithm. The user defined
   qinit function is called to compute value of conserved variables at
   the Gaussian quadrature points. These are then used to project the
   initial conditions on the basis functions

   Parameters
   ----------

   rd [in] - Initialized Run_Data object
   q [out] - Array of conserved variables [1-mbc:mx+mbc, 1:meqn, 1:sp_order]

*/
void
rkdg_initialize(Run_Data& rd, FArray<double>& q)
{
    double ec, wc, xcell;

    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;
    double xlower = rd.xlower;
    double dx = rd.dx;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work; // pointer to RKDG workspace
    int order = rd.sp_order; // order of Gaussian quarature
    
    // array to store cell centers 
    FArray<double> xloc( Range(1-mbc,mx+mbc), 0.0);
    // array to communicate with user's qinit
    FArray<double> qtemp( Range(1-mbc,mx+mbc, 1,meqn), 0.0);
    
    // loop over each abscissa calling qinit and accumulating the
    // result in q
    for(int c=1; c<=order; c++)
    {
        ec = rws->x(c); // eta of quadrature point
        wc = rws->w(c); // weight of quadrature point

        // compute cell coordinate corresponding to ec
        for(int i=1; i<=mx; i++)
        {
            xcell = xlower + (i-0.5)*dx; // cell center coordinate
            xloc(i) = ATE(ec,xcell,dx); // x-coordinate of ec
            
        }
        // get initial conditions at xloc
        rd.qinit(rd,xloc,qtemp);

        // accumulate result into q
        for(int i=1; i<=mx; i++)
            for(int m=1; m<=meqn; m++)
                for(int cc=1; cc<=order; cc++)
                    q(i,m,cc) += wc*rws->pmx(cc,c)*qtemp(i,m);

    }
    // Normalize the coefficients
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int cc=1; cc<=order; cc++)
                q(i,m,cc) = 0.5*(2.0*(cc-1.0)+1)*q(i,m,cc);
}
