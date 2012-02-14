#include "miniwarpx.h"
#include "rkdg_algo.h"
#include <math.h>

/**
   Computes the spatial part (right hand) of the system of hyperbolic
   equations
   
   dq/dt = -df/dx + s

   The 'rhs' is then used in a Runge-Kutta time stepping routine to
   advance the solution in time.
 */
static
void
rkdg_rhs(Run_Data& rd, FArray<double>& q, double tcurr, FArray<double>& rhs)
{
    int sgn;
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mbc = rd.mbc;
    int order = rd.sp_order;

    double xcell,xlower = rd.xlower;
    double ec,wc;
    double dx = rd.dx;
    double dx1 = 1./dx;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    rd.tcurrent = tcurr;

    //**
    // Step 1: Solve Reimann problem at cell interface
    //**

    // compute conserved variables and fluxes at left cell interfaces

    for(int i=1-mbc; i<=mx+mbc; i++)
        rd.xcoords(i) = xlower + i*dx; // left edge coordinates
    rkdg_eval_expansion_left_edge(rd,rws->ql,q);
    rd.flux(rd,rws->fl,rws->ql);

    // compute conserved variables and fluxes at right cell interfaces
        
    for(int i=1-mbc; i<=mx+mbc; i++)
        rd.xcoords(i) = xlower + (i+1.)*dx; // right edge coordinates
    rkdg_eval_expansion_right_edge(rd,rws->qr,q);
    rd.flux(rd,rws->fr,rws->qr);

    if(rd.edge_splitting == f_wave)
    { // we are using f-waves
        
        // compute the jump in flux and store it in df 
        for(int i=2-mbc; i<=mx+mbc; i++)
            for(int m=1; m<=meqn; m++)
                rws->df(i,m) = rws->fl(i,m) - rws->fr(i-1,m);
    }
    else
    { // we are using q-waves

        // compute the jump in conserved variables and store it in df
        for(int i=2-mbc; i<=mx+mbc; i++)
            for(int m=1; m<=meqn; m++)
                rws->df(i,m) = rws->ql(i,m) - rws->qr(i-1,m);
    }

    // solve Reimann problem at cell interface
    rd.rp(rd,rws->ql,rws->qr,rws->df,rws->wave,rws->s,rws->amdq,rws->apdq);

    // using solution of Reimann problem compute interface flux at
    // each edge
    for(int i=2-mbc; i<mx+mbc; i++)
        for(int m=1; m<=meqn; m++)
            rws->fedge(i,m) = 0.5*(rws->fl(i,m)+rws->fr(i-1,m)) -
                0.5*(rws->apdq(i,m)-rws->amdq(i,m));

    //**
    // Step 2: loop over each quadrature point, computing flux and
    // accumulating the result in rhs
    //**

    // initialize 'rhs' array: 'rhs' is initialized with the
    // difference of the edge fluxes computed from the Reimann solver
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
        {
            sgn = 1;
            for(int cc=1; cc<=order; cc++)
            {
                rhs(i,m,cc) = -dx1*(rws->fedge(i+1,m) - sgn*rws->fedge(i,m));
                sgn *= -1; // odd coefficients have -ve sign, while even have +ve
            }
        }

    FArray<double> qtemp = rws->df; // reuse array 'df' as it is no longer needed
    FArray<double> src = rws->fl; // reuse array 'fl' as it is no longer needed
    for(int c=1; c<=order; c++)
    {
        ec = rws->x(c); // eta of quadrature point
        wc = rws->w(c); // weight of quadrature point

        // compute cell coordinate corresponding to ec
        for(int i=1; i<=mx; i++)
        {
            xcell = xlower + (i-0.5)*dx; // cell center coordinate
            rd.xcoords(i) = ATE(ec,xcell,dx); // x-coordinate of ec
        }

        // compute conserved variable at quadrature point
        rkdg_eval_expansion(rd,qtemp,q,c);
        // compute flux at quadrature point
        rd.flux(rd,rws->f,qtemp);

        if (rd.has_source)
        { // source terms are present
 
            // compute source at quadrature point
            rd.src(rd,src,qtemp);
            // accumulate flux and src into 'rhs'
            for(int i=1; i<=mx; i++)
                for(int m=1; m<=meqn; m++)
                    for(int cc=1; cc<=order; cc++)
                        rhs(i,m,cc) += dx1*wc*rws->ppmx(cc,c)*rws->f(i,m)
                            + 0.5*wc*rws->pmx(cc,c)*src(i,m);
        }
        else
        { // no source terms are present

            // accumulate flux into 'rhs'
            for(int i=1; i<=mx; i++)
                for(int m=1; m<=meqn; m++)
                    for(int cc=1; cc<=order; cc++)
                        rhs(i,m,cc) += dx1*wc*rws->ppmx(cc,c)*rws->f(i,m);
        }
    }
    // Normalize the coefficients
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int cc=1; cc<=order; cc++)
                rhs(i,m,cc) = rhs(i,m,cc)/rws->Cconst(cc);

}

/**
   Runge-Kutta 1st order (forward Euler) scheme time advance
   algorithm. This is unstable for all spatial orders greater than 1
   and so should never be used, except for testing.
 */
static
void
rkdg_rk_order_1(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int ncoeffs = rd.ncoeffs;
    double dtdx = dt/rd.dx;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    //
    // Runge-Kutta Step 1
    //
    rkdg_rhs(rd,q,tcurr,rws->rhs); // compute RHS of equations
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int c=1; c<=ncoeffs; c++)
                q(i,m,c) = q(i,m,c) + dt*rws->rhs(i,m,c);

    // compute cfl number used
    cfla = 0.0;
    for(int i=1; i<=mx+1; i++)
        for(int m=1; m<=rd.mwave; m++)
             cfla = dmax(cfla, dtdx*fabs(rws->s(i,m)));
}

/**
   Runge-Kutta 2nd order TVD scheme time advance algorithm
 */
static
void
rkdg_rk_order_2_tvd(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int ncoeffs = rd.ncoeffs;
    double dtdx = dt/rd.dx;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    // apply limiters before RK step 1 if needed
    if(rd.dg_limiters != 0)
        rkdg_limiter(rd,q,dt);

    //
    // Runge-Kutta Step 1
    //
    rkdg_rhs(rd,q,tcurr,rws->rhs); // compute RHS of equations
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int c=1; c<=ncoeffs; c++)
                rws->q1(i,m,c) = q(i,m,c) + dt*rws->rhs(i,m,c);

    // compute cfl number used
    cfla = 0.0;
    for(int i=1; i<=mx+1; i++)
        for(int m=1; m<=rd.mwave; m++)
             cfla = dmax(cfla, dtdx*fabs(rws->s(i,m)));
    
    // check if we need to continue with computations
    if (cfla>rd.cflm)
        // no need to continue: although this step will be ultimately
        // redone it is best to return at this point as this avoids
        // problems with the RP being called again in the later RK stages
        return;

    // apply boundary conditions to get correct values in the ghost cells
    bc(rd,rws->q1);

    // apply limiters before RK step 2 if needed
    if(rd.dg_limiters != 0)
        rkdg_limiter(rd,rws->q1,dt);

    //
    // Runge-Kutta Step 2
    //
    rkdg_rhs(rd,rws->q1,tcurr,rws->rhs); // compute RHS of equation
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int c=1; c<=ncoeffs; c++)
                q(i,m,c) = 0.5*(q(i,m,c) + rws->q1(i,m,c) + dt*rws->rhs(i,m,c));

}

/**
   Runge-Kutta 3rd order TVD scheme time advance algorithm
 */
static
void
rkdg_rk_order_3_tvd(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int ncoeffs = rd.ncoeffs;
    double dtdx = dt/rd.dx;
    double c3b4 = 3./4.;
    double c2b3 = 2./3.;
    double c1b3 = 1./3.;
    double c1b4 = 1./4.;

    RKDG_Workspace *rws = (RKDG_Workspace*) rd.work;

    // apply limiters before RK step 1 if needed
    if(rd.dg_limiters != 0)
        rkdg_limiter(rd,q,dt);

    //
    // Runge-Kutta Step 1
    //
    rkdg_rhs(rd,q,tcurr,rws->rhs); // compute RHS of equations
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int c=1; c<=ncoeffs; c++)
                rws->q1(i,m,c) = q(i,m,c) + dt*rws->rhs(i,m,c);
    
    // compute cfl number used
    cfla = 0.0;
    for(int i=1; i<=mx+1; i++)
        for(int m=1; m<=rd.mwave; m++)
             cfla = dmax(cfla, dtdx*fabs(rws->s(i,m)));

    // check if we need to continue with computations
    if (cfla>rd.cflm)
        // no need to continue: although this step will be ultimately
        // redone it is best to return at this point as this avoids
        // problems with the RP being called again in the later RK stages
        return;

    // apply boundary conditions to get correct values in the ghost cells
    bc(rd,rws->q1);

    // apply limiters before RK step 2
    if(rd.dg_limiters != 0)
        rkdg_limiter(rd,rws->q1,dt);

    //
    // Runge-Kutta Step 2
    //
    rkdg_rhs(rd,rws->q1,tcurr,rws->rhs); // compute RHS of equation
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int c=1; c<=ncoeffs; c++)
                rws->q1(i,m,c) = c3b4*q(i,m,c) + c1b4*rws->q1(i,m,c)
                    + c1b4*dt*rws->rhs(i,m,c);

    // apply boundary conditions to get correct values in the ghost cells
    bc(rd,rws->q1);

    // apply limiters before RK step 3
    if(rd.dg_limiters != 0)
        rkdg_limiter(rd,rws->q1,dt);

    //
    // Runge-Kutta Step 3
    //
    rkdg_rhs(rd,rws->q1,tcurr,rws->rhs); // compute RHS of equation
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            for(int c=1; c<=ncoeffs; c++)
                q(i,m,c) = c1b3*q(i,m,c) + c2b3*rws->q1(i,m,c)
                    + c2b3*dt*rws->rhs(i,m,c);
}

/**
   Runge-Kutta 4nd order scheme time advance algorithm. This is NOT a
   TVD scheme
 */
static
void
rkdg_rk_order_4(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
}

/**
   Advances the solution of the hyperbolic conservation law by 'dt'
   using the Runge-Kutta Discontinous Galerkin method.

   Parameters
   ----------

   rd [in]    - Input data for simulation
   q [in/out] - On input, solution at start of step, on output,
                solution at end of step.
   tcurr [in] - Current time
   dt [in]    - Time step
   cfla [out] - CFL number used in this step.

 */
void 
rkdg_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
    if (rd.rk_order==1)
        rkdg_rk_order_1(rd,q,tcurr,dt,cfla);
    else if (rd.rk_order==2)
        rkdg_rk_order_2_tvd(rd,q,tcurr,dt,cfla);
    else if(rd.rk_order==3)
        rkdg_rk_order_3_tvd(rd,q,tcurr,dt,cfla);
    else if(rd.rk_order==4)
        rkdg_rk_order_4(rd,q,tcurr,dt,cfla);
    else
    { /* bad bad very bad */ }
}
