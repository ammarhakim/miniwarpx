#include "miniwarpx.h"
#include "euler.h"

#include <stdio.h>
#include <stdlib.h>
#include "dlapack_lite.h"
#include "blas_lite.h"

/*
  This function is called after each time step. Users can modify it to
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
after_step(Run_Data& rd, FArray<double>& q, double t, double dt)
{
    int INFO, LDA, LDB, M, N, NRHS;
    int *IPIV;
    char TRANS;
    int mx  = rd.mx;
    int mbc = rd.mbc;
    int meqn= rd.meqn;
    FArray<double> I(Range(1,5, 1,5),0.0);
    FArray<double> SM(Range(1,5, 1,5),0.0);
    FArray<double> A(Range(1,5, 1,5),0.0);
    FArray<double> B(Range(1,5, 1,1),0.0);
    // qbym
    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double qbym = ev->qbym;

    // generate identity matrix, I
    for (int i=1; i<=meqn; i++)
    {
        for (int j=1; j<=meqn; j++)
        {
            if (i==j) I(i,j)=1.0;
        }
    }

    for (int i=1-mbc; i<=mx+mbc; i++)
    {

        //generate source matrix, SM for source = SM*q
        SM(2,3) = qbym*ev->bf(i,3);
        SM(2,4) = -qbym*ev->bf(i,2);
        SM(3,2) = -qbym*ev->bf(i,3);
        SM(3,4) = qbym*ev->bf(i,1);
        SM(4,2) = qbym*ev->bf(i,2);
        SM(4,3) = -qbym*ev->bf(i,1);

        // define B as the q variables for the system A x = B
        for (int eqi=1; eqi<=meqn; eqi++)
        {
            B(eqi,1) = q(i,eqi,1);
            for (int eqj=1; eqj<=meqn; eqj++)
            {
                // define A as I-dt*SM
                A(eqi,eqj) = I(eqi,eqj) - dt*SM(eqi,eqj);
            }
        }
 
        M = N = LDA = LDB = A.dims(0);
        // allocate memory for permutation array
        IPIV = (int*) malloc( sizeof(int)*M );
        
        // compute LU factorization
        dgetrf_(&M, &N, A.data(), &LDA, IPIV, &INFO);

        // solve systems
        TRANS = 'N';
        NRHS = B.dims(1);
        // make call
        dgetrs_(&TRANS, &N, &NRHS, A.data(), &LDA, IPIV, B.data(), &LDB, &INFO);

        for (int eq=1; eq<=meqn; eq++)
        {
            q(i,eq,1) = B(eq,1);
        }
        free(IPIV);
    }
}
