#include "miniwarpx.h"
#include "wave_algo.h"
#include <math.h>
#include <iostream>

/**
   Advances the solution of the hyperbolic conservation law by 'dt'
   using the High Resolution Wave-Propagation method.

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
wave_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla)
{
    int mx = rd.mx;
    int meqn = rd.meqn;
    int mwave = rd.mwave;
    int mbc = rd.mbc;

    double xlower = rd.xlower;
    double dx = rd.dx;
    double dthalf = dt/2.0;
    double dtdxavg;

    FArray<double> xc(Range(1-mbc,mx+mbc),0.0);

    WAVE_Workspace *wws = (WAVE_Workspace*) rd.work;

    FArray<double> dtdx = wws->dtdx; // pointer to already allocated array

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
    // Source term advance: If we are using Strang splitting solve
    // dq/dt = s over half a time step
    //
    if (rd.has_source && (rd.source_splitting==2))
    {
        wave_source_advance(rd,q,tcurr,dthalf);
        // extend solution into boundary cells: is this really needed?
        bc(rd,q);
    }

    //
    // Solve homogeneous part of equation
    //

    //**
    // Step 1: Solve Reimann problem at cell interface
    //**

    for(int i=1-mbc; i<=mx+mbc; i++)
        for(int m=1; m<=meqn; m++)
            // simply copy the cell average into wws->q
            wws->q(i,m) = q(i,m,1);

    // compute the vector to be split into waves at each edge
    if(rd.edge_splitting==f_wave)
    {
        // we are using f-waves so compute jump in fluxes at each edge
        // compute conserved variables and fluxes at left cell interfaces

        rd.flux(rd,wws->fl,wws->q);

        // compute conserved variables and fluxes at right cell interfaces
        rd.flux(rd,wws->fr,wws->q);


        // compute the jump in flux and store it in fs
        for(int i=2-mbc; i<=mx+mbc; i++)
            for(int m=1; m<=meqn; m++)
                wws->fs(i,m) = wws->fl(i,m) - wws->fr(i-1,m);
    }
    else
    {
        // we are using q-waves so compute jump in q at each edge

        // compute the jump in conserved variables and store it in fs
        for(int i=2-mbc; i<=mx+mbc; i++)
            for(int m=1; m<=meqn; m++)
                wws->fs(i,m) = wws->q(i,m) - wws->q(i-1,m);
    }

    // solve Reimann problem at cell interface
    rd.rp(rd,wws->q,wws->q,wws->fs,wws->wave,wws->s,wws->amdq,wws->apdq);

    //**
    // Step 2: Use Reimann solution to compute first-order Gudunov updates
    //**

    for(int i=1; i<=mx+1; i++)
        for(int m=1; m<=meqn; m++)
        {
            q(i,m,1) += -dtdx(i)*wws->apdq(i,m);
            q(i-1,m,1) += - dtdx(i-1)*wws->amdq(i,m);
        }

    //**
    // Step 3: Compute CFL number for this step
    //**
    cfla = 0.0;
    for(int i=1; i<=mx; i++)
        for(int mw=1; mw<=mwave; mw++)
            cfla = dmax(cfla, dtdx(i)*wws->s(i,mw), -dtdx(i-1)*wws->s(i,mw));

    if (rd.wv_order==2)
    { // we are using second order Lax-Wendroff method

        //**
        // Step 4: Compute high resolution (second-order) corrections to fluxes
        //**
        
        // initialize correction fluxes
        for(int i=1-mbc; i<=mx+mbc; i++)
            for(int m=1; m<=meqn; m++)
                wws->fs(i,m) = 0.0;

        // apply limiters to waves
        wave_limiter(rd,wws->wave,wws->s);

        // compute second order corrections based on type of
        // edge_splitting being used
        if(rd.edge_splitting==f_wave)
        { // f-waves

            // compute second order correction fluxes
            for(int i=1; i<=mx+1; i++)
            {
                dtdxavg = 0.5*(dtdx(i-1)+dtdx(i)); // use "average" time step
                for(int m=1; m<=meqn; m++)
                    for(int mw=1; mw<=mwave; mw++)
                        wws->fs(i,m) += 0.5*dsign(1.0,wws->s(i,mw))
                            *(1.0-fabs(wws->s(i,mw))*dtdxavg)*wws->wave(i,m,mw);
            }
        }
        else
        { // q-waves

            // compute second order correction fluxes
            for(int i=1; i<=mx+1; i++)
            {
                dtdxavg = 0.5*(dtdx(i-1)+dtdx(i)); // use "average" time step
                for(int m=1; m<=meqn; m++)
                    for(int mw=1; mw<=mwave; mw++)
                        wws->fs(i,m) += 0.5*fabs(wws->s(i,mw))
                            *(1.0-fabs(wws->s(i,mw))*dtdxavg)*wws->wave(i,m,mw);
            }

        }

        // correct conserved variables
        for(int i=1; i<=mx; i++)
            for(int m=1; m<=meqn; m++)
                q(i,m,1) += -dtdx(i)*(wws->fs(i+1,m) - wws->fs(i,m));
    }
    // homoegenous part of equation solved


    //
    // Source term advance: If we are using Strang splitting solve
    // dq/dt = s over half a time step
    //
    if (rd.has_source && (rd.source_splitting==2))
        wave_source_advance(rd,q,tcurr+dt/2,dthalf);

    //
    // Source term advance: If we are using Gudonov splitting solve
    // dq/dt = s over full time step
    //
    if (rd.has_source && (rd.source_splitting==1))
        wave_source_advance(rd,q,tcurr,dt);
}
