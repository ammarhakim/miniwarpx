#include <iostream>
#include <fstream>
#include <math.h>

#include "farray.h"
#include "miniwarpx.h"

static
void
fa_copy(const FArray<double>& src, FArray<double>& dest)
{
    for(int i=src.start(0); i<=src.end(0); i++)
        for(int j=src.start(1); j<=src.end(1); j++)
            for(int k=src.start(2); k<=src.end(2); k++)
                dest(i,j,k) = src(i,j,k);
}

/**
   Advances the solution from time 'tcurr' to 'tend'. On entry 'q'
   should contain solution at 'tstart'. On exit 'q' contains solution
   at 'tend'.

   Parameters
   ----------

   rd [in/out] - Input data for this simulation
   q [in/out]  - On input must contain solution at time 'tcurr'. On
                 output contains solution at time 'tend'.
   tcurr [in]  - Time at start of step
   tend [in]   - Time at end of step

 */
void 
advance(Run_Data& rd, FArray<double>& q, double tcurr, double tend)
{

    double told, cfla, dtlast;
    FArray<double> qcopy(q.range(), 0.0);
    
    double t = tcurr;
    double dt = rd.nv[0];
    double cfl = rd.cfl;
    double cflm = rd.cflm;
    double nstep = 1;

    // loop, advancing solution using adaptive time steps
    while(1)
    {
        told = t;

        dtlast = dt;

        // adjust dt to hit tend exactly if we are near the end of the
        // computation
        if (told+dt>tend) dt = tend-told;

        rd.tcurrent = t; // set current time before bc is called
        // apply boundary condition to get correct values in ghost cells
        bc(rd,q);

        // copy solution q into qcopy in case we have to repeat this
        // step
        fa_copy(q,qcopy);
        
      redo:
        t = told + dt; // time at end of step

        // call before_step to give user chance to do things before
        // the solution is advanced
        before_step(rd,q,told);

        // advance solution by dt
        rd.step(rd,q,told,dt,cfla);

        // call after_step to give user chance to do things after the
        // solution is advanced
        after_step(rd,q,t,dt);

        if (rd.verbose)
            std::cout << "MINWARPX Step " << nstep 
                      << " Time " << t 
                      << " dt " << dt 
                      << " CFL " << cfla 
                      << std::endl;

        // check if Courant number was too large
        if(cfla>cflm)
        {
            // courant number was too large: go back and redo this step

            dt = dt*cfl/cfla; // new time step
            t = told; // reset starting time
            fa_copy(qcopy,q); // copy solution at previous step back into q
            std::cout << "** MINIWARPX rejecting step... Courant number " << cfla 
                      <<  " too large" 
                      << std::endl;
            
            goto redo;
        }
        // adjust time step to ensure we are getting desired CFL number
        dt = dt*cfl/cfla;

        nstep += 1;
        if(t >= tend) break; // break if we are done
    }
    rd.nv[0] = dtlast; // copy last time step
}
