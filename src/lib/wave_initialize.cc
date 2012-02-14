#include "miniwarpx.h"

/**
   Set initial conditions for the Wave Propagation algorithm. The user
   defined qinit function is called to compute value of conserved
   variables at the cell centers.

   Parameters
   ----------

   rd [in] - Initialized Run_Data object
   q [out] - Array of conserved variables [1-mbc:mx+mbc, 1:meqn, 1]

*/
void
wave_initialize(Run_Data& rd, FArray<double>& q)
{

    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;

    // array to store cell centers 
    FArray<double> xloc_c( Range(1-mbc,mx+mbc), 0.0);
    FArray<double> xloc_p( Range(1-mbc,mx+mbc), 0.0);
    // array to communicate with user's qinit
    FArray<double> qtemp( Range(1-mbc,mx+mbc, 1,meqn), 0.0);

    double xlower = rd.xlower;
    double dx = rd.dx;

    // compute cell center coordinates in computational space
    for(int i=1; i<=mx; i++)
        xloc_c(i) = xlower + (i-0.5)*dx;
    // transform computational space coordinates into physical space
    grid_transform(rd,xloc_c,xloc_p);

    // set initial conditions by calling user's qinit
    rd.qinit(rd,xloc_p,qtemp);
    for(int i=1; i<=mx; i++)
        for(int m=1; m<=meqn; m++)
            q(i,m,1) = qtemp(i,m);

    if(rd.has_kappa)
        // set capacity function
        rd.setkappa(rd,xloc_p,rd.kappa);
}
