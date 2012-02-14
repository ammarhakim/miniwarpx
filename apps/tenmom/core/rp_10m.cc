#include <stdio.h>
#include <iostream>
#include <math.h>
#include "miniwarpx.h"

extern
void
eigensystem(const Run_Data& rd, double *ql, double *qr, double *ev, double **lev, double **rev);

void 
rp_10m(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
       FArray<double>& wave, FArray<double>& s, 
       FArray<double>& amdq, FArray<double>& apdq)
{
    // Reimann solver using full eigensystem and arithmetic averages

    int mx = rd.mx;
    int mbc = rd.mbc;

    double qrt[10], qlt[10], ev[10], **rev, **lev;
    double a[10];

    rev = new double*[10];
    lev = new double*[10];
    for(int i=0; i<10; i++)
    {
        rev[i] = new double[10];
        lev[i] = new double[10];
    }

    // compute waves
    for(int i=2-mbc; i<=mx+mbc; i++)
    {
        // compute eigenvalues and right/left eigenvectors using averaged conserved variables
        for(int m=1; m<=10; m++)
        {
            qlt[m-1] = qr(i-1,m);
            qrt[m-1] = ql(i,m);
        }
        eigensystem(rd,qlt,qrt,ev,lev,rev);

        // project df onto left eigenvectors
        for(int mw=0; mw<10; mw++)
        {
            a[mw] = 0.0;
            for(int m=0; m<10; m++)
                a[mw] += lev[mw][m]*df(i,m+1);
        }
        // compute waves and wave speeds
        for(int mw=0; mw<10; mw++)
        {
            for(int m=0; m<10; m++)
                wave(i,m+1,mw+1) = a[mw]*rev[m][mw];
            s(i,mw+1) = ev[mw]; // wave speed
        }
    }
    // compute fluctuations
    eval_fluctuations(rd,wave,s,amdq,apdq);

    for(int i=0; i<10; i++)
    {
        delete [] rev[i];
        delete [] lev[i];
    }
    delete [] rev;
    delete [] lev;
}
