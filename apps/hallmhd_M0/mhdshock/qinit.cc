#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "hallmhd.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;
    double nl_e,pl_e;
    double nl_i,ul_i,vl_i,wl_i,pl_i;
    double nr_e,pr_e;
    double nr_i,ur_i,vr_i,wr_i,pr_i;
    double El_e,El_i;
    double Er_e,Er_i;
    double bxr,byr,bzr,exr,eyr,ezr;
    double bxl,byl,bzl,exl,eyl,ezl;

    Hallmhd_Vars *hv = (Hallmhd_Vars*) rd.mvar;

    double gas_gamma  = hv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = hv->mi;
    double qi = hv->qi;
    double c0 = hv->c0;

    // set left initial state
    // for ions
    nl_i = 1.0;
    ul_i = 0.0;
    vl_i = 0.0;
    wl_i = 0.0;
    pl_i = 0.5e-4;
    // for electrons
    nl_e = 1.0;
    pl_e = 0.5e-4;
    // for EM field
    bxl = 0.75e-2;
    byl = 0.0;
    bzl = -1.0e-2;
    exl = 0.0;
    eyl = 0.0;
    ezl = 0.0;

    hv->ur    = 0.1;
//sqrt(pl_i/(mi*nl_i));
    hv->betar = 1.0;
//pl_i/(bxl*bxl);
    hv->rli   = 0.01;
//mi/(qi*bxl*L)*hv->ur;

    // set right initial state...
    // for ions
    nr_i = 0.125;
    ur_i = 0.0;
    vr_i = 0.0;
    wr_i = 0.0;
    pr_i = 0.05e-4;
    // for electrons
    nr_e = 0.125;
    pr_e = 0.05e-4;
    // for EM field
    bxr = 0.75e-2;
    byr = 0.0;
    bzr = 1.0e-2;
    exr = 0.0;
    eyr = 0.0;
    ezr = 0.0;

    // calculate energy...
    // for left state electrons and ions
    El_i = pl_i/gas_gamma1 + 0.5*mi*nl_i*(pow(ul_i,2)+pow(vl_i,2)+pow(wl_i,2));
    El_e = pl_e/gas_gamma1;
    // for right state electrons and ions
    Er_i = pr_i/gas_gamma1 + 0.5*mi*nr_i*(pow(ur_i,2)+pow(vr_i,2)+pow(wr_i,2));
    Er_e = pr_e/gas_gamma1;

    double ur = hv->ur;
    double sloc = 0.5;
    double sloc1 = 0.6;
    double dx = sloc1-sloc;
    for(int i=1-rd.mbc; i<=rd.mx+rd.mbc; i++)
    {
        xcell = xloc(i);
        if(xcell < sloc)
        {
            // left state
            q(i,1)  = nl_i;
            q(i,2)  = nl_i*ul_i;
            q(i,3)  = nl_i*vl_i;
            q(i,4)  = nl_i*wl_i;
            q(i,5)  = El_i;
            q(i,6)  = nl_e;
            q(i,7)  = El_e;
            q(i,8)  = (ur*ur)/(c0*c0)*exl;
            q(i,9)  = (ur*ur)/(c0*c0)*eyl;
            q(i,10) = (ur*ur)/(c0*c0)*ezl;
            q(i,11) = bxl;
            q(i,12) = byl;
            q(i,13) = bzl;
        }
        else if(xcell >= sloc && xcell<=sloc1)
        {
            // for maccormack treat like line with slope
            q(i,1)  = (nr_i-nl_i)/dx*xcell + nl_i-(nr_i-nl_i)/dx*sloc;
            q(i,2)  = ((nr_i-nl_i)/dx*xcell + nl_i-(nr_i-nl_i)/dx*sloc)*
                ((ur_i-ul_i)/dx*xcell + ul_i-(ur_i-ul_i)/dx*sloc);
            q(i,3)  = ((nr_i-nl_i)/dx*xcell + nl_i-(nr_i-nl_i)/dx*sloc)*
                ((vr_i-vl_i)/dx*xcell + vl_i-(vr_i-vl_i)/dx*sloc);
            q(i,4)  = ((nr_i-nl_i)/dx*xcell + nl_i-(nr_i-nl_i)/dx*sloc)*
                ((wr_i-wl_i)/dx*xcell + wl_i-(wr_i-wl_i)/dx*sloc);
            q(i,5)  = (Er_i-El_i)/dx*xcell + El_i-(Er_i-El_i)/dx*sloc;
            q(i,6)  = (nr_e-nl_e)/dx*xcell + nl_e-(nr_e-nl_e)/dx*sloc;
            q(i,7)  = (Er_e-El_e)/dx*xcell + El_e-(Er_e-El_e)/dx*sloc;
            q(i,8)  = (ur*ur)/(c0*c0)*((exr-exl)/dx*xcell + exl-(exr-exl)/dx*sloc);
            q(i,9)  = (ur*ur)/(c0*c0)*((eyr-eyl)/dx*xcell + eyl-(eyr-eyl)/dx*sloc);
            q(i,10) = (ur*ur)/(c0*c0)*((ezr-ezl)/dx*xcell + ezl-(ezr-ezl)/dx*sloc);
            q(i,11) = (bxr-bxl)/dx*xcell + bxl-(bxr-bxl)/dx*sloc;
            q(i,12) = (byr-byl)/dx*xcell + byl-(byr-byl)/dx*sloc;
            q(i,13) = (bzr-bzl)/dx*xcell + bzl-(bzr-bzl)/dx*sloc;
        }
        else
        {
            // right state
            q(i,1)  = nr_i;
            q(i,2)  = nr_i*ur_i;
            q(i,3)  = nr_i*vr_i;
            q(i,4)  = nr_i*wr_i;
            q(i,5)  = Er_i;
            q(i,6)  = nr_e;
            q(i,7)  = Er_e;
            q(i,8)  = (ur*ur)/(c0*c0)*exr;
            q(i,9)  = (ur*ur)/(c0*c0)*eyr;
            q(i,10) = (ur*ur)/(c0*c0)*ezr;
            q(i,11) = bxr;
            q(i,12) = byr;
            q(i,13) = bzr;
        }
    }
}
