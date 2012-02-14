#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "tenmom.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    Tenmom_Vars *tfv = (Tenmom_Vars*) rd.mvar;

    double mi = tfv->mi;
    double me = tfv->me;

    double ratio_me_mi = me/mi;
    double sloc,xcell;
    double rhol_e,ul_e,vl_e,wl_e,pl_e;
    double rhol_i,ul_i,vl_i,wl_i,pl_i;
    double rhor_e,ur_e,vr_e,wr_e,pr_e;
    double rhor_i,ur_i,vr_i,wr_i,pr_i;
    double exl,eyl,ezl,bxl,byl,bzl;
    double exr,eyr,ezr,bxr,byr,bzr;

    sloc = 0.5; // location of shock

    // set left initial state...
    // for electrons
    rhol_e = 1.0*ratio_me_mi;
    ul_e = 0.e0;
    vl_e = 0.e0;
    wl_e = 0.e0;
    pl_e = 0.5e-4;
    //  for ions
    rhol_i = 1.0;
    ul_i = 0.e0;
    vl_i = 0.e0;
    wl_i = 0.e0;
    pl_i = 0.5e-4;
    //  for EM field
    bxl = 0.75e-2;
    byl = 0.0;
    bzl = -1.0e-2;
    exl = 0.e0;
    eyl = 0.e0;
    ezl = 0.e0;

    //  set right initial state...
    // for electrons
    rhor_e = 0.125*ratio_me_mi;
    ur_e = 0.e0;
    vr_e = 0.e0;
    wr_e = 0.e0;
    pr_e = 0.05e-4;
    // for ions
    rhor_i = 0.125;
    ur_i = 0.e0;
    vr_i = 0.e0;
    wr_i = 0.e0;
    pr_i = 0.05e-4;
    // for EM field
    bxr = 0.75e-2;
    byr = 0.0;
    bzr = 1.0e-2;
    exr = 0.e0;
    eyr = 0.e0;
    ezr = 0.e0;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        if(xcell < sloc)
        {
            // left state
            q(i,1) = rhol_e; // electrons
            q(i,2) = rhol_e*ul_e;
            q(i,3) = rhol_e*vl_e;
            q(i,4) = rhol_e*wl_e;
            q(i,5) = pl_e + rhol_e*ul_e*ul_e;
            q(i,6) = rhol_e*ul_e*vl_e;
            q(i,7) = rhol_e*ul_e*wl_e;
            q(i,8) = pl_e + rhol_e*vl_e*vl_e;
            q(i,9) = rhol_e*vl_e*wl_e;
            q(i,10) = pl_e + rhol_e*wl_e*wl_e;

            q(i,11) = rhol_i; // ions
            q(i,12) = rhol_i*ul_i;
            q(i,13) = rhol_i*vl_i;
            q(i,14) = rhol_i*wl_i;
            q(i,15) = pl_i + rhol_i*ul_i*ul_i;
            q(i,16) = rhol_i*ul_i*vl_i;
            q(i,17) = rhol_i*ul_i*wl_i;
            q(i,18) = pl_i + rhol_i*vl_i*vl_i;
            q(i,19) = rhol_i*vl_i*wl_i;
            q(i,20) = pl_i + rhol_i*wl_i*wl_i;

            q(i,21) = exl; // EM fields 
            q(i,22) = eyl;
            q(i,23) = ezl;
            q(i,24) = bxl;
            q(i,25) = byl;
            q(i,26) = bzl;
        }
        else
        {
            // right state
            q(i,1) = rhor_e; // electrons
            q(i,2) = rhor_e*ur_e;
            q(i,3) = rhor_e*vr_e;
            q(i,4) = rhor_e*wr_e;
            q(i,5) = pr_e + rhor_e*ur_e*ur_e;
            q(i,6) = rhor_e*ur_e*vr_e;
            q(i,7) = rhor_e*ur_e*wr_e;
            q(i,8) = pr_e + rhor_e*vr_e*vr_e;
            q(i,9) = rhor_e*vr_e*wr_e;
            q(i,10) = pr_e + rhor_e*wr_e*wr_e;

            q(i,11) = rhor_i; // ions
            q(i,12) = rhor_i*ur_i;
            q(i,13) = rhor_i*vr_i;
            q(i,14) = rhor_i*wr_i;
            q(i,15) = pr_i + rhor_i*ur_i*ur_i;
            q(i,16) = rhor_i*ur_i*vr_i;
            q(i,17) = rhor_i*ur_i*wr_i;
            q(i,18) = pr_i + rhor_i*vr_i*vr_i;
            q(i,19) = rhor_i*vr_i*wr_i;
            q(i,20) = pr_i + rhor_i*wr_i*wr_i;

            q(i,21) = exr; // EM fields
            q(i,22) = eyr;
            q(i,23) = ezr;
            q(i,24) = bxr;
            q(i,25) = byr;
            q(i,26) = bzr;
        }
    }
}
