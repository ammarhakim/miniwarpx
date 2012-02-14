#include "tenmom.h"
#include "miniwarpx.h"
#include "LAPACK/dlapack_lite.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

void
eigensystem(const Run_Data& rd, double *ql, double *qr, double *ev, double **levc, double **revc)
{
    double rl,ul,vl,wl,pxxl,pxyl,pxzl,pyyl,pyzl,pzzl;
    double rr,ur,vr,wr,pxxr,pxyr,pxzr,pyyr,pyzr,pzzr;

    // compute primitive variables at left and right states
    rl = ql[0];
    ul = ql[1]/rl;
    vl = ql[2]/rl;
    wl = ql[3]/rl;
    pxxl = ql[4]-rl*ul*ul;
    pxyl = ql[5]-rl*ul*vl;
    pxzl = ql[6]-rl*ul*wl;
    pyyl = ql[7]-rl*vl*vl;
    pyzl = ql[8]-rl*vl*wl;
    pzzl = ql[9]-rl*wl*wl;

    rr = qr[0];
    ur = qr[1]/rr;
    vr = qr[2]/rr;
    wr = qr[3]/rr;
    pxxr = qr[4]-rr*ur*ur;
    pxyr = qr[5]-rr*ur*vr;
    pxzr = qr[6]-rr*ur*wr;
    pyyr = qr[7]-rr*vr*vr;
    pyzr = qr[8]-rr*vr*wr;
    pzzr = qr[9]-rr*wr*wr;

    // average the left and right states
    double r,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz;
    
    r = 0.5*(rl+rr);
    u = 0.5*(ul+ur);
    v = 0.5*(vl+vr);
    w = 0.5*(wl+wr);
    pxx = 0.5*(pxxr+pxxl);
    pxy = 0.5*(pxyr+pxyl);
    pxz = 0.5*(pxzr+pxzl);
    pyy = 0.5*(pyyr+pyyl);
    pyz = 0.5*(pyzr+pyzl);
    pzz = 0.5*(pzzr+pzzl);

    double cxx2,cxy2,cxz2,cyy2,cyz2,czz2;
    cxx2 = pxx/r;
    cxy2 = pxy/r;
    cxz2 = pxz/r;
    cyy2 = pyy/r;
    cyz2 = pyz/r;
    czz2 = pzz/r;

    // compute eigenvalues: arranged in increasing order
    double cxx = sqrt(cxx2);
    double r3 = sqrt(3.0);

    ev[0] = u-r3*cxx;
    ev[1] = ev[2] = u-cxx;
    ev[3] = ev[4] = ev[5] = ev[6] = u;
    ev[7] = ev[8] = u+cxx;
    ev[9] = u+r3*cxx;

    double lev[10][10], rev[10][10];

    // initialize eigenvector matrices as many are zero
    for (int i=0; i<10; ++i)
        for (int j=0; j<10; ++j)
        {
            lev[i][j] = 0.0;
            rev[i][j] = 0.0;
        }

    // compute right eigenvectors of primitive variables: these are
    // stored in columns

    // r1
    rev[0][0] = 1.0;
    rev[1][0] = -r3*cxx/r;
    rev[2][0] = -r3*cxy2/(r*cxx);
    rev[3][0] = -r3*cxz2/(r*cxx);
    rev[4][0] = 3.0*cxx2;
    rev[5][0] = 3.0*cxy2;
    rev[6][0] = 3.0*cxz2;
    rev[7][0] = (cxx2*cyy2+2*cxy2*cxy2)/cxx2;
    rev[8][0] = (cxx2*cyz2+2*cxy2*cxz2)/cxx2;
    rev[9][0] = (cxx2*czz2+2*cxz2*cxz2)/cxx2;

    // r2
    rev[2][1] = -cxx/r;
    rev[5][1] = cxx2;
    rev[7][1] = 2*cxy2;
    rev[8][1] = cxz2;

    // r3
    rev[3][2] = -cxx/r;
    rev[6][2] = cxx2;
    rev[8][2] = cxy2;
    rev[9][2] = 2*cxz2;

    // r4
    rev[0][3] = 1.0;

    // r5
    rev[7][4] = cxx2;

    // r6
    rev[8][5] = cxx2;

    // r7
    rev[9][6] = cxx2;

    // r8
    rev[3][7] = cxx/r;
    rev[6][7] = cxx2;
    rev[8][7] = cxy2;
    rev[9][7] = 2*cxz2;

    // r9
    rev[2][8] = cxx/r;
    rev[5][8] = cxx2;
    rev[7][8] = 2*cxy2;
    rev[8][8] = cxz2;

    // r10
    rev[0][9] = 1.0;
    rev[1][9] = r3*cxx/r;
    rev[2][9] = r3*cxy2/(r*cxx);
    rev[3][9] = r3*cxz2/(r*cxx);
    rev[4][9] = 3.0*cxx2;
    rev[5][9] = 3.0*cxy2;
    rev[6][9] = 3.0*cxz2;
    rev[7][9] = (cxx2*cyy2+2*cxy2*cxy2)/cxx2;
    rev[8][9] = (cxx2*cyz2+2*cxy2*cxz2)/cxx2;
    rev[9][9] = (cxx2*czz2+2*cxz2*cxz2)/cxx2;

    // compute left eigenvectors of primitive variables: these are
    // stored as row vectors

    // l1
    lev[0][1] = -r3*r/(6*cxx);
    lev[0][4] = 1.0/(6*cxx2);

    // l2
    lev[1][1] = r*cxy2/(2*cxx2*cxx);
    lev[1][2] = -r/(2*cxx);
    lev[1][4] = -cxy2/(2*cxx2*cxx2);
    lev[1][5] = 1.0/(2*cxx2);

    // l3
    lev[2][1] = r*cxz2/(2*cxx2*cxx);
    lev[2][3] = -r/(2*cxx);
    lev[2][4] = -cxz2/(2*cxx2*cxx2);
    lev[2][6] = 1.0/(2*cxx2);

    // l4
    lev[3][0] = 1.0;
    lev[3][4] = -1.0/(3*cxx2);

    // l5
    lev[4][4] = (4*cxy2*cxy2-cxx2*cyy2)/(3*cxx2*cxx2*cxx2);
    lev[4][5] = -2*cxy2/(cxx2*cxx2);
    lev[4][7] = 1/cxx2;

    // l6
    lev[5][4] = (4*cxy2*cxz2-cxx2*cyz2)/(3*cxx2*cxx2*cxx2);
    lev[5][5] = -cxz2/(cxx2*cxx2);
    lev[5][6] = -cxy2/(cxx2*cxx2);
    lev[5][8] = 1/cxx2;

    // l7
    lev[6][4] = (4*cxz2*cxz2-cxx2*czz2)/(3*cxx2*cxx2*cxx2);
    lev[6][6] = -2*cxz2/(cxx2*cxx2);
    lev[6][9] = 1/cxx2;

    // l8
    lev[7][1] = -r*cxz2/(2*cxx2*cxx);
    lev[7][3] = r/(2*cxx);
    lev[7][4] = -cxz2/(2*cxx2*cxx2);
    lev[7][6] = 1/(2*cxx2);

    // l9
    lev[8][1] = -r*cxy2/(2*cxx2*cxx);
    lev[8][2] = r/(2*cxx);
    lev[8][4] = -cxy2/(2*cxx2*cxx2);
    lev[8][5] = 1/(2*cxx2);

    // l10
    lev[9][1] = r3*r/(6*cxx);
    lev[9][4] = 1/(6*cxx2);
    
    // now construct the M and Minv matrices to compute the conserved
    // variable eigensystem
    FArray<double> M(Range(0,9, 0,9)), Minv(Range(0,9, 0,9));

    // initialize matrices as most entries are just zero
    for (int i=0; i<10; ++i)
        for (int j=0; j<10; ++j)
        {
            M(i,j) = 0.0;
            Minv(i,j) = 0.0;
        }

    // initialize M

    // row 1
    M(0,0) = 1.0;

    // row 2
    M(1,0) = u;
    M(1,1) = r;

    // row 3
    M(2,0) = v;
    M(2,2) = r;

    // row 4
    M(3,0) = w;
    M(3,3) = r;

    // row 5
    M(4,0) = u*u;
    M(4,1) = 2*r*u;
    M(4,4) = 1.0;
    
    // row 6
    M(5,0) = u*v;
    M(5,1) = r*v;
    M(5,2) = r*u;
    M(5,5) = 1.0;

    // row 7
    M(6,0) = u*w;
    M(6,1) = r*w;
    M(6,3) = r*u;
    M(6,6) = 1.0;

    // row 8
    M(7,0) = v*v;
    M(7,2) = 2*r*v;
    M(7,7) = 1.0;

    // row 9
    M(8,0) = v*w;
    M(8,2) = r*w;
    M(8,3) = r*v;
    M(8,8) = 1.0;

    // row 10
    M(9,0) = w*w;
    M(9,3) = 2*r*w;
    M(9,9) = 1.0;

    // copy M into Minv
    for (int i=0; i<10; ++i)
        for (int j=0; j<10; ++j)
            Minv(i,j) = M(i,j);

    // compute minv using LAPACK rouinte dtrtri_
    int info, neq = 10;
    dtrtri_("L", "N", &neq, Minv.data(), &neq, &info);
    // check if inversion worked
    if (info != 0)
    {
        std::cout << "Inversion of lower-triangular matrix failed" << std::endl;
        exit(1);
    }

    // compute conserved variable eigensystem by doing matrix multiply
    for (int i=0; i<10; ++i)
        for (int j=0; j<10; ++j)
        {
            revc[i][j] = 0.0;
            levc[i][j] = 0.0;
        }

    for (int i=0; i<10; ++i)
        for (int j=0; j<10; ++j)
            for (int k=0; k<10; ++k)
            {
                revc[i][j] += M(i,k)*rev[k][j];
                levc[i][j] += lev[i][k]*Minv(k,j);
            }
}
