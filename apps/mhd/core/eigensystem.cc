/**
 * Computes Ideal MHD equation eigensystem.
 */

#include <math.h>
#include "miniwarpx.h"
#include "mhd.h"

/**
   Computes the eigensystem of the ideal-MHD equations. 

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

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double g2g1 = (gas_gamma-2)/(gas_gamma-1);

    double rhol,ul,vl,wl,prl,bxl,byl,bzl;
    double rhor,ur,vr,wr,prr,bxr,byr,bzr;

    rhol = ql[0];    rhor = qr[0];
    ul = ql[1]/rhol; ur = qr[1]/rhor;
    vl = ql[2]/rhol; vr = qr[2]/rhor;
    wl = ql[3]/rhol; wr = qr[3]/rhor;

    bxl = ql[5];     bxr = qr[5];
    byl = ql[6];     byr = qr[6];
    bzl = ql[7];     bzr = qr[7];

    prl = gas_gamma1*(ql[4] - 0.5*rhol*(ul*ul+vl*vl+wl*wl) - 0.5*(bxl*bxl+byl*byl+bzl*bzl));
    prr = gas_gamma1*(qr[4] - 0.5*rhor*(ur*ur+vr*vr+wr*wr) - 0.5*(bxr*bxr+byr*byr+bzr*bzr));

    // compute averaged primitive variables
    double rho = 0.5*(rhor+rhol);
    double u = 0.5*(ul+ur);
    double v = 0.5*(vl+vr);
    double w = 0.5*(wl+wr);

    double bx = 0.5*(bxl+bxr);
    double by = 0.5*(byl+byr);
    double bz = 0.5*(bzl+bzr);

    double pr = 0.5*(prl+prr);

    double btot2 = bx*bx + by*by + bz*bz;

    // compute speeds
    double a = sqrt(gas_gamma*pr/rho); // speed of sound
    double ca = sqrt(bx*bx/rho); // Alfven speed

    double t1 = a*a + btot2/rho;
    double t2 = sqrt(t1*t1 - 4*a*a*bx*bx/rho);

    double cf = sqrt(0.5*(t1 + t2)); // Fast magnetosonic speed
    double cs = sqrt(0.5*(t1 - t2)); // Slow magnetosonic speed

    // compute eigenvalues: these are arranged from slowest to fastest
    ev[0] = u-cf;
    ev[1] = u-ca;
    ev[2] = u-cs;
    ev[3] = u;
    ev[4] = u+cs;
    ev[5] = u+ca;
    ev[6] = u+cf;

    
    // define various quantities to make computing eigenvectors easier
    double b1, b2, b3, th1, th2, alf, als;
    double mufm, mufp, musm, musp;
    double t3,t4;

    if (bx==0) 
        b1=0;
    else
        b1 = (bx>0) ? 1 : -1;
    if ((by==bz) && (by==0.0))
        b2 = b3 = 1.0/sqrt(2.0);
    else
    {
        b2 = by/sqrt(by*by+bz*bz);
        b3 = bz/sqrt(by*by+bz*bz);
    }

    if ((by==bz) && (by==0.0) && (ca*ca == a*a))
        alf = als = 1.0;
    else
    {
        alf = sqrt((cf*cf-ca*ca)/(cf*cf-cs*cs));
        als = sqrt((cf*cf-a*a)/(cf*cf-cs*cs));
    }

    t1 = alf*cf*cf/gas_gamma1;
    t2 = alf*cf*u;
    t3 = als*ca*b1*(b2*v+b3*w);
    t4 = g2g1*alf*(cf*cf-a*a);

    mufm = t1 - t2 + t3 + t4;
    mufp = t1 + t2 - t3 + t4;

    t1 = als*cs*cs/gas_gamma1;
    t2 = als*cs*u;
    t3 = alf*a*b1*(b2*v+b3*w);
    t4 = g2g1*als*(cs*cs-a*a);

    musm = t1 - t2 - t3 + t4;
    musp = t1 + t2 + t3 + t4;

    th1 = alf*alf*a*a*(cf*cf-g2g1*a*a) + als*als*cf*cf*(cs*cs-g2g1*a*a);
    th1 = 0.5/th1;
    th2 = alf*alf*cf*a*b1 + als*als*cs*ca*b1;
    th2 = 0.5/th2;

    double u2 = u*u + v*v + w*w;
    // compute right eigenvectors: these are stored as column vectors

    // R1 
    rev[0][0] = alf;
    rev[1][0] = alf*(u-cf);
    rev[2][0] = alf*v+als*b1*b2*ca;
    rev[3][0] = alf*w+als*b1*b3*ca;
    rev[4][0] = 0.5*alf*u2 + mufm;
    rev[5][0] = 0.0;
    rev[6][0] = als*b2*cf/sqrt(rho);
    rev[7][0] = als*b3*cf/sqrt(rho);

    // R7: this is similar to R1 so can make more efficient
    rev[0][6] = alf;
    rev[1][6] = alf*(u+cf);
    rev[2][6] = alf*v-als*b1*b2*ca;
    rev[3][6] = alf*w-als*b1*b3*ca;
    rev[4][6] = 0.5*alf*u2 + mufp;
    rev[5][6] = 0.0;
    rev[6][6] = als*b2*cf/sqrt(rho);
    rev[7][6] = als*b3*cf/sqrt(rho);

    // R2
    rev[0][1] = 0.0;
    rev[1][1] = 0.0;
    rev[2][1] = b1*b3;
    rev[3][1] = -b1*b2;
    rev[4][1] = (b3*v-b2*w)*b1;
    rev[5][1] = 0.0;
    rev[6][1] = b3/sqrt(rho);
    rev[7][1] = -b2/sqrt(rho);

    // R6: this is similar to R2 so can make more efficient
    rev[0][5] = 0.0;
    rev[1][5] = 0.0;
    rev[2][5] = -b1*b3;
    rev[3][5] = b1*b2;
    rev[4][5] = -(b3*v-b2*w)*b1;
    rev[5][5] = 0.0;
    rev[6][5] = b3/sqrt(rho);
    rev[7][5] = -b2/sqrt(rho);

    // R4
    rev[0][3] = 1.0;
    rev[1][3] = u;
    rev[2][3] = v;
    rev[3][3] = w;
    rev[4][3] = 0.5*u2;
    rev[5][3] = 0.0;
    rev[6][3] = 0.0;
    rev[7][3] = 0.0;

    // R3
    rev[0][2] = als;
    rev[1][2] = als*(u-cs);
    rev[2][2] = als*v-alf*b1*b2*a;
    rev[3][2] = als*w-alf*b1*b3*a;
    rev[4][2] = 0.5*als*u2 + musm;
    rev[5][2] = 0.0;
    rev[6][2] = -alf*b2*a*a/(cf*sqrt(rho));
    rev[7][2] = -alf*b3*a*a/(cf*sqrt(rho));

    // R5
    rev[0][4] = als;
    rev[1][4] = als*(u+cs);
    rev[2][4] = als*v+alf*b1*b2*a;
    rev[3][4] = als*w+alf*b1*b3*a;
    rev[4][4] = 0.5*als*u2 + musp;
    rev[5][4] = 0.0;
    rev[6][4] = -alf*b2*a*a/(cf*sqrt(rho));
    rev[7][4] = -alf*b3*a*a/(cf*sqrt(rho));

    // compute left eigenvectors: these are stored as row vectors

    // L1
    lev[0][0] = 0.5*th1*alf*a*a*u2 + th2*(alf*a*u*b1 - als*cs*(b2*v+b3*w));
    lev[0][1] = -th1*alf*a*a*u - th2*alf*a*b1;
    lev[0][2] = -th1*alf*a*a*v + th2*als*cs*b2;
    lev[0][3] = -th1*alf*a*a*w + th2*als*cs*b3;
    lev[0][4] = th1*alf*a*a;
    lev[0][5] = 0.0;
    lev[0][6] = th1*sqrt(rho)*als*b2*cf*(cs*cs-g2g1*a*a);
    lev[0][7] = th1*sqrt(rho)*als*b3*cf*(cs*cs-g2g1*a*a);

    // L7
    lev[6][0] = 0.5*th1*alf*a*a*u2 - th2*(alf*a*u*b1 - als*cs*(b2*v+b3*w));
    lev[6][1] = -th1*alf*a*a*u + th2*alf*a*b1;
    lev[6][2] = -th1*alf*a*a*v - th2*als*cs*b2;
    lev[6][3] = -th1*alf*a*a*w - th2*als*cs*b3;
    lev[6][4] = th1*alf*a*a;
    lev[6][5] = 0.0;
    lev[6][6] = th1*sqrt(rho)*als*b2*cf*(cs*cs-g2g1*a*a);
    lev[6][7] = th1*sqrt(rho)*als*b3*cf*(cs*cs-g2g1*a*a);

    // L2
    lev[1][0] = -0.5*b1*(b3*v-b2*w);
    lev[1][1] = 0.0;
    lev[1][2] = 0.5*b1*b3;
    lev[1][3] = -0.5*b1*b2;
    lev[1][4] = 0.0;
    lev[1][5] = 0.0;
    lev[1][6] = 0.5*sqrt(rho)*b3;
    lev[1][7] = -0.5*sqrt(rho)*b2;

    // L6
    lev[5][0] = 0.5*b1*(b3*v-b2*w);
    lev[5][1] = 0.0;
    lev[5][2] = -0.5*b1*b3;
    lev[5][3] = 0.5*b1*b2;
    lev[5][4] = 0.0;
    lev[5][5] = 0.0;
    lev[5][6] = 0.5*sqrt(rho)*b3;
    lev[5][7] = -0.5*sqrt(rho)*b2;

    // L4
    lev[3][0] = 1.0 - th1*(alf*alf*a*a + als*als*cf*cf)*u2;
    lev[3][1] = 2.0*th1*(alf*alf*a*a + als*als*cf*cf)*u;
    lev[3][2] = 2.0*th1*(alf*alf*a*a + als*als*cf*cf)*v;
    lev[3][3] = 2.0*th1*(alf*alf*a*a + als*als*cf*cf)*w;
    lev[3][4] = -2.0*th1*(alf*alf*a*a + als*als*cf*cf);
    lev[3][5] = 0.0;
    lev[3][6] = 2.0*th1*sqrt(rho)*alf*als*b2*cf*(cf*cf-cs*cs);
    lev[3][7] = 2.0*th1*sqrt(rho)*alf*als*b3*cf*(cf*cf-cs*cs);

    // L3
    lev[2][0] = 0.5*th1*als*cf*cf*u2 + th2*(als*ca*u*b1 + alf*cf*(b2*v+b3*w));
    lev[2][1] = -th1*als*cf*cf*u - th2*als*ca*b1;
    lev[2][2] = -th1*als*cf*cf*v - th2*alf*cf*b2;
    lev[2][3] = -th1*als*cf*cf*w - th2*alf*cf*b3;
    lev[2][4] = th1*als*cf*cf;
    lev[2][5] = 0.0;
    lev[2][6] = -th1*sqrt(rho)*alf*b2*cf*(cf*cf-g2g1*a*a);
    lev[2][7] = -th1*sqrt(rho)*alf*b3*cf*(cf*cf-g2g1*a*a);

    // L5
    lev[4][0] = 0.5*th1*als*cf*cf*u2 - th2*(als*ca*u*b1 + alf*cf*(b2*v+b3*w));
    lev[4][1] = -th1*als*cf*cf*u + th2*als*ca*b1;
    lev[4][2] = -th1*als*cf*cf*v + th2*alf*cf*b2;
    lev[4][3] = -th1*als*cf*cf*w + th2*alf*cf*b3;
    lev[4][4] = th1*als*cf*cf;
    lev[4][5] = 0.0;
    lev[4][6] = -th1*sqrt(rho)*alf*b2*cf*(cf*cf-g2g1*a*a);
    lev[4][7] = -th1*sqrt(rho)*alf*b3*cf*(cf*cf-g2g1*a*a);

}
