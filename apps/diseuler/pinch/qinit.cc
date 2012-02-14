#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "miniwarpx.h"
#include "euler.h"

double
prf(double r)
{
    return -512./3.*pow(r,6) + 12.*pow(r,4) - 1./4.*pow(r,2);
}

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    double pr;

    Euler_Vars *ev = (Euler_Vars*) rd.mvar;
    double gas_gamma  = ev->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    double lambda = ev->qbym;

    double a = 1./8.;
    // compute parameters for pressure profile
    double alpha = 0.1;
    double pa = prf(a);
    double p0 = lambda*pa/(alpha-1.);

    for(int i=1; i<=rd.mx; i++)
    {
        // radial coordinate
        double r = xloc(i);

        if(r<a)
        {
            pr = p0 + lambda*prf(r);

            q(i,1) = 1.0;
            q(i,2) = 0.0; 
            q(i,3) = 0.0;
            q(i,4) = q(i,1)*(1-64.*pow(r,2));
            q(i,5) = pr/gas_gamma1 + 0.5*pow(q(i,4),2)/q(i,1);
        }
        else
        {
            pr = p0 + lambda*prf(a);

            q(i,1) = 1.0;
            q(i,2) = 0.0; 
            q(i,3) = 0.0;
            q(i,4) = 0.0;
            q(i,5) = pr/gas_gamma1;
        }
    }

    for(int i=1; i<=rd.mx+rd.mbc; i++)
    {
        // radial coordinate
        double r = xloc(i);

        if(r<a)
        {
            // initialize magnetic field
            ev->bf(i,1) = 0.0;
            ev->bf(i,2) = 0.5*r - 16.*pow(r,3);
            ev->bf(i,3) = 0.0;
        }
        else
        {
            // initialize magnetic field
            ev->bf(i,1) = 0.0;
            ev->bf(i,2) = a*(0.5*a - 16.*pow(a,3))/r;
            ev->bf(i,3) = 0.0;
        }
    }

    // apply axis boundary condition to B field
          
    for(int ibc=1; ibc<=rd.mbc; ibc++)
    {
        // radial terms are flipped in sign
        ev->bf(1-ibc,1) = -a*ev->bf(ibc,1);
        // phi terms are flipped in sign
        ev->bf(1-ibc,2) = -a*ev->bf(ibc,2);
        // z terms are copied across axis
        ev->bf(1-ibc,3) = a*ev->bf(ibc,3);
    }

    std::ofstream fout("b");
    for(int i=1; i<=rd.mx; i++)
    {
        for(int c=1; c<=3; c++)
            fout << ev->bf(i,c) << " ";
        fout << std::endl;
    }
}
