#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double rhol_e,ul_e,vl_e,wl_e,pl_e;
    double rhol_i,ul_i,vl_i,wl_i,pl_i;
    double rhor_e,ur_e,vr_e,wr_e,pr_e;
    double rhor_i,ur_i,vr_i,wr_i,pr_i;
    double El_e,El_i;
    double Er_e,Er_i;
    double bxr,byr,bzr,exr,eyr,ezr;
    double bxl,byl,bzl,exl,eyl,ezl;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = tfv->mi;
    double me = tfv->me;

    double ratio_me_mi = me/mi;
    // set left initial state
    // for electrons
    rhol_e = 1.0*ratio_me_mi;
    ul_e = 0.0;
    vl_e = 0.0;
    wl_e = 0.0;
    pl_e = 0.5e-4;
    // for ions
    rhol_i = 1.0;
    ul_i = 0.0;
    vl_i = 0.0;
    wl_i = 0.0;
    pl_i = 0.5e-4;
    // for EM field
    bxl = 0.75e-2;
    byl = 0.0;
    bzl = -1.0e-2;
    exl = 0.0;
    eyl = 0.0;
    ezl = 0.0;

    // set right initial state...
    // for electrons
    rhor_e = 0.125*ratio_me_mi;
    ur_e = 0.0;
    vr_e = 0.0;
    wr_e = 0.0;
    pr_e = 0.05e-4;
    // for ions
    rhor_i = 0.125;
    ur_i = 0.0;
    vr_i = 0.0;
    wr_i = 0.0;
    pr_i = 0.05e-4;
    // for EM field
    bxr = 0.75e-2;
    byr = 0.0;
    bzr = 1.0e-2;
    exr = 0.0;
    eyr = 0.0;
    ezr = 0.0;

    // calculate energy...
    // for left state electrons and ions
    El_e = pl_e/gas_gamma1 + 0.5*rhol_e*(pow(ul_e,2)+pow(vl_e,2)+pow(wl_e,2));
    El_i = pl_i/gas_gamma1 + 0.5*rhol_i*(pow(ul_i,2)+pow(vl_i,2)+pow(wl_i,2));
    // for right state electrons and ions
    Er_e = pr_e/gas_gamma1 + 0.5*rhor_e*(pow(ur_e,2)+pow(vr_e,2)+pow(wr_e,2));
    Er_i = pr_i/gas_gamma1 + 0.5*rhor_i*(pow(ur_i,2)+pow(vr_i,2)+pow(wr_i,2));

    double sloc = 0.5;
    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);
        if(xcell < sloc)
        {
            // left state
            q(i,1)  = rhol_e;
            q(i,2)  = rhol_e*ul_e;
            q(i,3)  = rhol_e*vl_e;
            q(i,4)  = rhol_e*wl_e;
            q(i,5)  = El_e;
            q(i,6)  = rhol_i;
            q(i,7)  = rhol_i*ul_i;
            q(i,8)  = rhol_i*vl_i;
            q(i,9)  = rhol_i*wl_i;
            q(i,10) = El_i;
            q(i,11) = exl;
            q(i,12) = eyl;
            q(i,13) = ezl;
            q(i,14) = bxl;
            q(i,15) = byl;
            q(i,16) = bzl;
        }
        else
        {
            // right state
            q(i,1)  = rhor_e;
            q(i,2)  = rhor_e*ur_e;
            q(i,3)  = rhor_e*vr_e;
            q(i,4)  = rhor_e*wr_e;
            q(i,5)  = Er_e;
            q(i,6)  = rhor_i;
            q(i,7)  = rhor_i*ur_i;
            q(i,8)  = rhor_i*vr_i;
            q(i,9)  = rhor_i*wr_i;
            q(i,10) = Er_i;
            q(i,11) = exr;
            q(i,12) = eyr;
            q(i,13) = ezr;
            q(i,14) = bxr;
            q(i,15) = byr;
            q(i,16) = bzr;
        }
    }
}
