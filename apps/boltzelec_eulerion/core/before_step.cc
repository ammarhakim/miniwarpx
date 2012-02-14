#include "miniwarpx.h"
#include "euler.h"

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include "dlapack_lite.h"
#include "blas_lite.h"

/*
  This function is called before each time step. Users can modify it to
  do anything which needs to be done after the step is taken.

  Solves systems of linear equations A*x = B

  Inputs
  ------
  A - Coefficient matrix
  B - Righthand sides stored as columns of B

  Output
  ------

  B - Each column of B is replaced by the solution vector for the
  corresponding righthand side
  
*/

void
before_step(Run_Data& rd, FArray<double>& q, double t)
{
    int INFO, LDA, LDB, M, N, NRHS;
    char TRANS;
    int mx    = rd.mx;
    int mbc   = rd.mbc;
    double dx = rd.dx;
    double n0 = 1.0;
    Euler_Vars *ev = (Euler_Vars*) rd.mvar;

    FArray<double> B(Range(1,mx, 1,1),0.0);
    FArray<double> niedge(Range(1,mx),0.0);

    // compute edge values for ni
    for (int i=1; i<=mx; i++)
        niedge(i) = (q(i-1,1,1)+q(i,1,1))/2.0;

    // solving for phi iteratively 
    double tol = 1.e-10;
    double err = 1.0;
    int num = 1;
    int iter = 10;
    
    while (err>tol && num<iter)
    {

        // apply periodic boundary conditions
        ev->phi(0) = ev->phi(1) - (ev->phi(mx) - ev->phi(mx-1));
        ev->phi(mx+1) = ev->phi(mx) + (ev->phi(2) - ev->phi(1));
        cout << "Possion residual " << err << endl;

        // compute the RHS of Laplace's eqn to solve iteratively for phi
        for (int i=2; i<mx; i++)
        {
            ev->dF(i,i)   = -2./(dx*dx)-n0*exp(-ev->phi(i));
            ev->dF(i,i-1) = 1./(dx*dx);
            ev->dF(i,i+1) = 1./(dx*dx);
        }
        ev->dF(1,1)       = -2./(dx*dx)-n0*exp(-ev->phi(1));
        ev->dF(1,2)       = 1./(dx*dx);
        ev->dF(mx,mx) = -2./(dx*dx)-n0*exp(-ev->phi(mx));
        ev->dF(mx,mx-1)   = 1./(dx*dx);

        for (int i=1; i<=mx; i++)
            B(i,1) = 1./(dx*dx)*(ev->phi(i-1) - 2.0*ev->phi(i) + ev->phi(i+1))
                - (niedge(i) - n0*exp(-ev->phi(i)));
        
        M = N = LDA = LDB = ev->dF.dims(0);        
        // compute LU factorization
        dgetrf_(&M, &N, ev->dF.data(), &LDA, ev->IPIV.data(), &INFO);

        // solve systems
        TRANS = 'N';
        NRHS = B.dims(1);
        // make call
        dgetrs_(&TRANS, &N, &NRHS, ev->dF.data(), &LDA, ev->IPIV.data(), B.data(), &LDB, &INFO);

        // compute 2 norm to check for convergence
        err=0.;
        for (int j=1; j<=mx; j++)
            err += pow(B(j,1),2);
        err = sqrt(err);
        num = num+1;

        for (int i=1; i<=mx; i++)
            ev->phi(i) = ev->phi(i) - B(i,1);

    }

    // compute the electric field (cell center values) using E=-grad(phi)
    for (int i=1; i<=mx; i++)
    {
        ev->ef(i,1) = (ev->phi(i+1)-ev->phi(i))/dx;
        ev->ef(i,2) = 0;
        ev->ef(i,3) = 0;
    }
}
