#include "miniwarpx.h"
#include "wave_algo.h"
#include "euler.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dlapack_lite.h"
#include "blas_lite.h"

/*
  Solves the ODE dq/dt = s from time 'tcurr' to 'tcurr+dt'. This
  particular routine uses a semi-implicit source term advance to solve the
  ODE. 

  Parameters
  ----------

  rd [in]     - Simulation parameters
  q  [in/out] - On input contains solution at time 'tcurr'. On output contains 
                solution at time 'tcurr+dt'.
  t  [in]     - Current time at which source term is being calculated
  dt [in]     - Time step to advace ODE by

 */
void 
wave_source_advance(Run_Data& rd, FArray<double>& q,double t, double dt)
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
    FArray<double> source(Range(1,5, 1,1),0.0);

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;

    // qbym
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
            source(eqi,1) = 0.0; // initialize source term
            for (int eqj=1; eqj<=meqn; eqj++)
            {
                // define A as I-dt/2*SM
                A(eqi,eqj) = I(eqi,eqj) - dt/2*SM(eqi,eqj);
                source(eqi,1) += (I(eqi,eqj)+dt/2*SM(eqi,eqj))*q(i,eqj,1);
            }
        }
 
        M = N = LDA = LDB = A.dims(0);
        // allocate memory for permutation array
        IPIV = (int*) malloc( sizeof(int)*M );
        
        // compute LU factorization
        dgetrf_(&M, &N, A.data(), &LDA, IPIV, &INFO);

        // solve systems
        TRANS = 'N';
        NRHS = source.dims(1);
        // make call
        dgetrs_(&TRANS, &N, &NRHS, A.data(), &LDA, IPIV, source.data(), &LDB, &INFO);

         for (int eqi=1; eqi<=meqn; eqi++)
         {
             //q(i,eqi,1) = q(i,eqi,1)+dt*source(eqi,1);
             q(i,eqi,1) = source(eqi,1);
         }

        free(IPIV);
    }
}
