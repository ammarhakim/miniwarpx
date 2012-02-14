#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "twofluid.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double xcell;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;

    double gas_gamma  = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double mi = tfv->mi;
    double me = tfv->me;
    double epsilon0 = tfv->epsilon0;
    double qe = tfv->qe;

    double pi = 3.141592654;

    double n0 = 1.0; // number density of electrons and ions
    double pr = 1.e-2; // small value of pressure: cold fluids
    double er;
    er = pr/gas_gamma1; // all energy is due to pressure

    for(int i=1; i<=rd.mx; i++)
    {
        // x-coordinate
        xcell = xloc(i);

        // uniform cold plasma with a EM pulse launched from left wall

        // electrons
        q(i,1) = me*n0;
        q(i,2) = 0.0;
        q(i,3) = 0.0;
        q(i,4) = 0.0;
        q(i,5) = er;

        // ions
        q(i,6)  = mi*n0;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        q(i,10) = er;

        // electric field
        q(i,11) = 0.0;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        // magnetic field
        q(i,14) = 0.0;
        q(i,15) = 0.0;
        q(i,16) = 0.0;
    }

}
