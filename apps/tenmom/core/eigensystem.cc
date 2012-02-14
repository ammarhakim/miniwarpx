/**
 * Computes Tenmoment equation eigensystem.
 */

#include <math.h>
#include "miniwarpx.h"

extern "C"
{
    void tm_sys_eig__(int *, double*, double*, double*, double*, double*);
}

/**
   Computes the eigensystem of the Tenmoment equations. 

   Parameters
   ----------

   rd - Run time information
   q - Conserved variables array
   ev - Eigenvalues
   lev - left eigenvectors as row vectors
   rev - right eigenvectors as column vectors
 */
void
eigensystem(const Run_Data& rd, double *ql, double *qr, double *ev, double **lev, double **rev)
{
    WxIndexer idx();
    FArray<double> lev_f(WxRange(0,9, 0,9)), rev_f(WxRange(0,9, 0,9));

    // call Fortran routine to compute eigensystem
    int meqn = 10;
    tm_sys_eig__(&meqn, ql, qr, ev, lev_f.data(), rev_f.data());
    // copy into appropriate locations
    for (int i=0; i<10; ++i)
        for (int j=0; j<10; ++j)
        {
            lev[i][j] = lev_f(i,j);
            rev[i][j] = rev_f(i,j);
        }
}
