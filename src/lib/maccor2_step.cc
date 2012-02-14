#include "miniwarpx.h"
#include "maccor_algo.h"
#include <math.h>
#include <iostream>

/**
   Advances the solution of the hyperbolic conservation law by 'dt'
   using the MacCormick 2nd order method.

   Parameters
   ----------

   rd [in]    - Input data for simulation
   q [in/out] - On input, solution at start of step, on output,
                solution at end of step.
   tcurr [in] - Current time
   dt [in]    - Time step
   cfla [out] - CFL number used in this step.

 */

static int nstep = 0;

void 
maccor2_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mbc = rd.mbc;

    double xlower = rd.xlower;
    double dx = rd.dx;

    FArray<double> xc(Range(1-mbc,mx+mbc),0.0);

    MACCOR2_Workspace *mws = (MACCOR2_Workspace*) rd.work;

    FArray<double> dtdx = mws->dtdx; // pointer to already allocated array
    double dthalf = dt/2.0;

    // set current time and cell-center coordinates in case flux and
    // src routines need them
    rd.tcurrent = tcurr; // set current time
    for(int i=1-mbc; i<=mx+mbc; i++)
        xc(i) = xlower + (i-0.5)*dx;
    // transform xc to physical domain
    grid_transform(rd,xc,rd.xcoords);

    // set dtdx array to take into account the capacity fucntion
    if(rd.has_kappa)
        // system has capacity function
        for(int i=1-mbc; i<=mx+mbc; i++)
            dtdx(i) = dt/(dx*rd.kappa(i));
    else
        // system does not have capacity function
        for(int i=1-mbc; i<=mx+mbc; i++)
            dtdx(i) = dt/dx;

    //
    // Source term advance: We are using Strang splitting so solve
    // dq/dt = s over half a time step
    //
    if (rd.has_source)
    {
        maccor2_source_advance(rd,q,tcurr,dthalf);
        // extend solution into boundary cells: is this really needed?
        bc(rd,q);
    }

    //
    // Solve homogeneous part of equation
    //

    for(int i=1-mbc; i<=mx+mbc; i++)
        for(int m=1; m<=meqn; m++)
            // simply copy the cell average into mws->q
            mws->q(i,m) = q(i,m,1);

    // two step 2nd order MacCormack scheme
    switch(nstep)
    {
        case 0:

            //
            // STEP 1: Predict solution 
            //
            
            for(int i=1-mbc; i<=mx+mbc; i++)
                for(int m=1; m<=meqn; m++)
                    // simply copy the cell average into mws->q
                    mws->q(i,m) = q(i,m,1);

            // compute flux at each grid point
            rd.flux(rd,mws->f,mws->q);
                
            // loop over all interior cells, predicting solution using
            // forward-differences
            for(int i=1; i<=mx; i++)
                for(int m=1; m<=meqn; m++)
                    mws->qp(i,m,1) = q(i,m,1) - dtdx(i)*(mws->f(i+1,m)-mws->f(i,m));
            
            // apply boundary conditions
            bc(rd,mws->qp);
            
            //
            // STEP 2: Correct solution
            //

            for(int i=1-mbc; i<=mx+mbc; i++)
                for(int m=1; m<=meqn; m++)
                    // simply copy the cell average into mws->q
                    mws->q(i,m) = mws->qp(i,m,1);
            
            // compute flux at each grid point
            rd.flux(rd,mws->f,mws->q);
            // compute maxmium wave speed at each grid point
            rd.maxs(rd,mws->s,mws->q);
            
            // loop over all interior cells, correcting solution using
            // backward-differences
            for(int i=1; i<=mx; i++)
                for(int m=1; m<=meqn; m++)
                    q(i,m,1) = 0.5*(q(i,m,1) + mws->qp(i,m,1) 
                                    - dtdx(i)*(mws->f(i,m)-mws->f(i-1,m)));

            // apply boundary conditions
            bc(rd,q);
            
            break;

        case 1:

            //
            // STEP 1: Predict solution 
            //
            
            for(int i=1-mbc; i<=mx+mbc; i++)
                for(int m=1; m<=meqn; m++)
                    // simply copy the cell average into mws->q
                    mws->q(i,m) = q(i,m,1);

            // compute flux at each grid point
            rd.flux(rd,mws->f,mws->q);
                
            // loop over all interior cells, predicting solution using
            // backward-differences
            for(int i=1; i<=mx; i++)
                for(int m=1; m<=meqn; m++)
                    mws->qp(i,m,1) = q(i,m,1) - dtdx(i)*(mws->f(i,m)-mws->f(i-1,m));
            
            // apply boundary conditions
            bc(rd,mws->qp);
            
            //
            // STEP 2: Correct solution
            //

            for(int i=1-mbc; i<=mx+mbc; i++)
                for(int m=1; m<=meqn; m++)
                    // simply copy the cell average into mws->q
                    mws->q(i,m) = mws->qp(i,m,1);
            
            // compute flux at each grid point
            rd.flux(rd,mws->f,mws->q);
            // compute maxmium wave speed at each grid point
            rd.maxs(rd,mws->s,mws->q);
            
            // loop over all interior cells, correcting solution using
            // forward-differences
            for(int i=1; i<=mx; i++)
                for(int m=1; m<=meqn; m++)
                    q(i,m,1) = 0.5*(q(i,m,1) + mws->qp(i,m,1) 
                                    - dtdx(i)*(mws->f(i+1,m)-mws->f(i,m)));

            // apply boundary conditions
            bc(rd,q);
            
            break;

    }

    // compute CFL number used
    cfla = 0.0;
    for(int i=1; i<=mx+1; i++) cfla = dmax(cfla,dtdx(i)*fabs(mws->s(i)));
    // homoegenous part of equation solved


    //
    // Source term advance: We are using Strang splitting solve so
    // dq/dt = s over half a time step
    //
    if (rd.has_source)
        maccor2_source_advance(rd,q,tcurr+dt/2,dthalf);

    nstep = (nstep+1) % 2; // toggle so we can switch between
                           // forward-backward and backward-forward
                           // stepping
}
