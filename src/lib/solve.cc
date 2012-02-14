#include <iostream>
#include <fstream>
#include <math.h>

#include "farray.h"
#include "miniwarpx.h"

/**
   Solves a system of hyperbolic conservation laws

   dq/dt + df/dx = s

   in one dimension. The algorithm used is either the High-Resolution
   Wave Propagation method or the Runge-Kutta Discontinuous Galerkin
   method.

   Parameters
   ----------

   rd [in/out] - Input data for simulation. This must be initialized
   correctly before a call to this function.

*/
void
solve(Run_Data& rd)
{
    // set run parameters from rd
    int meqn = rd.meqn;
    int mbc = rd.mbc;
    int mx = rd.mx;
    double tstart = rd.tstart, tend = rd.tend;
    int nout = rd.nout;
    rd.nv[0] = rd.dt;
    rd.nv[1] = rd.cfl;
    rd.nv[2] = rd.cflm;

    int failed;
    failed = 0;

    rd.dx = (rd.xupper-rd.xlower)/mx; // grid spacing

    // set algorithm specific data
    if (rd.algo == WAVE)
        // High-resolution wave propagation algorithm
        wave_setup(rd);
    else if (rd.algo == RKDG)
        // Runge-Kutta Discontinous Galerkin algorithm
        rkdg_setup(rd);
    else if(rd.algo == MACCOR2)
        // 2nd order MacCormick algorithm
        maccor2_setup(rd);
    else
    {/* bad bad very bad */ }

    // q: array of conserved variables
    FArray<double> q(Range(1-mbc,mx+mbc, 1,meqn, 1,rd.ncoeffs),0.0);

    // allocate memory for cell coordinates array
    rd.xcoords = FArray<double>(Range(1-mbc,mx+mbc), 0.0);

    // allocate memory for capacity function
    rd.kappa = FArray<double>(Range(1-mbc,mx+mbc), 1.0);

    // compute cell center coordinates in computational space
    FArray<double> xc( Range(1-mbc,mx+mbc), 0.0 );
    for(int i=1-mbc; i<=mx+mbc; i++)
        xc(i) = rd.xlower + (i-0.5)*rd.dx;
    // transform xc to physical domain
    grid_transform(rd,xc,rd.xcoords);
    // write out the coordinates to file frame.x 
    write_grid(rd,rd.xcoords);

    //**
    // Step 1: Set initial conditions and capacity function
    //**
    if (rd.algo==WAVE)
        // call initialization routine for WAVE
        wave_initialize(rd,q);
    else if(rd.algo == RKDG)
        // call initilization routine for RKDG
        rkdg_initialize(rd,q);
    else if(rd.algo == MACCOR2)
        // call initialization routine for MACCOR2
        maccor2_initialize(rd,q);
    else
    {/* bad bad very bad */ }

    //**
    // Step 2: Apply boundary conditions
    //**
    bc(rd,q);

    // write initial conditions to file
    out(rd,"frame.q0",q);
    std::cout << "MINIWARPX Frame " << 0 << " written to file" << std::endl;

    //**
    // Step 3: Advance solution in time
    //**
    std::cout << "Starting simulation..." << std::endl << std::endl;
    double tcurr = tstart; // starting time
    double tsize = (tend-tstart)/nout; // time between file output
    int frame = 1;
    for(int i=0; i<rd.nout; i++)
    {
        // advance solution by tsize
        advance(rd,q,tcurr,tcurr+tsize);

        // construct file name to write to
        char buff[100];
        sprintf(buff, "frame.q%d", frame);

        // write solution to file
        out(rd,buff,q);

        std::cout << "MINIWARPX Frame " << frame 
                  << " at time " << tcurr+tsize 
                  <<  " written to file" 
                  << std::endl;
        std::cout << std::endl;

        // advance tcurr and frame number
        tcurr += tsize;
        frame += 1;
    }
}
