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
    double epsilon0 = tfv->epsilon0;
    double mi = tfv->mi;
    double me = tfv->me;
    double qi = tfv->qi;
    double c0 = tfv->c0;

    double ratio_me_mi = me/mi;
    double rhoe,p,pe,pi,by,ex,mez;

    double n0    = 1.0;
    double B0    = 1.0;
    double deby  = sqrt(epsilon0*mi*c0*c0/(n0*qi*qi));
    double ni    = 1.0;
    double len   = 1.0;
    double tau   = len/c0;
    double rad   = 1.0;
    double Z     = 1.0;
    double wc    = qi*B0/mi;
    double alpha = 4*Z*Z*len*len*ni*ni/(deby*deby);
    double beta  = 4*Z*Z*wc*wc*tau*tau*ni*ni;
    double j0    = 1.0; //0.1*n0*qi*c0;
    double p0    = rad*rad*j0*j0;
    double nback = 100.0;
    double pback = 1.0;
//    double A, b1,b2,b3,b4,b5, c1,c2,c3,c4,c5, d;

    for(int i=1; i<=rd.mx; i++)
    {
        xcell = xloc(i);

//        if (xcell <= rad )
//        {

            rhoe = Z*ni*mi + deby*deby*p0/(2.*Z*ni*len*len*rad*rad)
                *(1.-3.*xcell*xcell/(rad*rad)+3.*pow(xcell,4)/(2.*pow(rad,4)));

            p  = p0*(5./48. - xcell*xcell/(4.*rad*rad) 
                       + 3.*pow(xcell,4)/(16.*pow(rad,4)) - pow(xcell,6)/(24.*pow(rad,6)));

            by = sqrt(p0/(4.*beta))*sqrt((p0/pow(rad,2)+alpha)*pow(xcell,2)/pow(rad,2)
                                    - (3.*p0/pow(rad,2)+alpha)*pow(xcell,4)/pow(rad,4)
                                    + 1./4.*(13.*p0/pow(rad,2)+alpha)*pow(xcell,6)/pow(rad,6)
                                    - 3./2.*p0/pow(rad,2)*pow(xcell,8)/pow(rad,8)
                                    + 1./4.*p0/pow(rad,2)*pow(xcell,10)/pow(rad,10));

            mez= -(wc*tau*deby*deby)/(len*len)
                * ((p0/(8.*beta*rad*by) * (2.*(p0/pow(rad,2)+alpha)*xcell/rad
                                      - 4.*(3.*p0/pow(rad,2)+alpha)*pow(xcell,3)/pow(rad,3)
                                      + 3./2.*(13.*p0/pow(rad,2)+alpha)*pow(xcell,5)/pow(rad,5)
                                      - 12.*p0/pow(rad,2)*pow(xcell,7)/pow(rad,7)
                                      + 5./2.*p0/pow(rad,2)*pow(xcell,9)/pow(rad,9)))
                   + by/xcell);

            ex = -p0/(2*Z*wc*tau*ni*rad)*(xcell/(2*rad)-3*pow(xcell,3)/(4*pow(rad,3))
                                          +pow(xcell,5)/(4*pow(rad,5)));

//         }
//         else
//         {
//             ne = nback + Z*ni + deby*deby*p0/(2.*Z*ni*len*len*rad*rad)*(1.-3.+3./2.);
//             p  = pback;

//             A  = p0/(8.*beta*rad);
//             b1 = 2.*(p0/(rad*rad)+alpha);
//             b2 = 4.*(3.*p0/(rad*rad)+alpha);
//             b3 = 3./2.*(13.*p0/(rad*rad)+alpha);
//             b4 = 12.*p0/(rad*rad);
//             b5 = 5./2.*p0/(rad*rad);

//             d  = sqrt(p0/(4*beta));
//             c1 = p0/(rad*rad)+alpha;
//             c2 = 3.*p0/(rad*rad)+alpha;
//             c3 = 1./4.*(13.*p0/(rad*rad)+alpha);
//             c4 = 3./2.*p0/(rad*rad);
//             c5 = 1./4.*p0/(rad*rad);

//             by = ((1./6.*(A*b1*rad+d*d*c1)/(pow(rad,2)*d*sqrt(c1/pow(rad,2))))*pow(rad,3)+(-1./40.*(2.*A*b2*c1*rad-A*c2*b1*rad+d*d*c1*c2)/(pow(rad,4)*d*sqrt(c1/pow(rad,2))*c1))*pow(rad,5)+(1./336.*(8.*A*b3*pow(c1,2)*rad-4.*A*c2*b2*c1*rad-4.*A*b1*c3*c1*rad+3.*A*b1*pow(c2,2)*rad+4.*d*d*pow(c1,2)*c3-d*d*c1*pow(c2,2))/(pow(rad,6)*d*sqrt(c1/pow(rad,2))*pow(c1,2)))*pow(rad,7)+(-1./1152.*(8.*d*d*pow(c1,3)*c4-4.*d*d*pow(c1,2)*c2*c3+d*d*c1*pow(c2,3)+16.*A*b4*pow(c1,3)*rad-8.*A*c2*b3*pow(c1,2)*rad-8.*A*b2*pow(c1,2)*c3*rad+6.*A*b2*c1*pow(c2,2)*rad-8.*A*b1*c4*pow(c1,2)*rad+12.*A*b1*c2*c3*c1*rad-5.*A*b1*pow(c2,3)*rad)/(pow(rad,8)*d*sqrt(c1/pow(rad,2))*pow(c1,3)))*pow(rad,9)+(1./14080.*(128.*A*b5*pow(c1,4)*rad-64.*A*c2*b4*pow(c1,3)*rad-64.*A*b3*pow(c1,3)*c3*rad+48.*A*b3*pow(c1,2)*pow(c2,2)*rad-64.*A*b2*pow(c1,3)*c4*rad+96.*A*b2*pow(c1,2)*c2*c3*rad-40.*A*b2*c1*pow(c2,3)*rad-120.*A*b1*pow(c2,2)*c3*c1*rad-64.*A*b1*c5*pow(c1,3)*rad+96.*A*b1*c2*c4*pow(c1,2)*rad+48.*A*b1*pow(c3,2)*pow(c1,2)*rad+35.*A*b1*pow(c2,4)*rad+64.*d*d*pow(c1,4)*c5-32.*d*d*pow(c1,3)*c2*c4-16.*d*d*pow(c1,3)*pow(c3,2)+24.*d*d*pow(c1,2)*pow(c2,2)*c3-5.*d*d*c1*pow(c2,4))/(pow(rad,10)*d*sqrt(c1/pow(rad,2))*pow(c1,4)))*pow(rad,11))/xcell + 0.40767;


//             mez= 0.0;
//             ex = 0.0;
//         }

        pe = 0.5*p;
        pi = pe;

        // initial conditions for electron fluid, ion fluid and e-m terms
        q(i,1)  = rhoe;
        q(i,2)  = 0.0;
        q(i,3)  = 0.0;
        q(i,4)  = mez;
        q(i,5)  = pe/gas_gamma1 + 0.5*(pow(q(i,2),2)+pow(q(i,3),2)
                                      +pow(q(i,4),2))/q(i,1);
        q(i,6)  = mi*ni;
        q(i,7)  = 0.0;
        q(i,8)  = 0.0;
        q(i,9)  = 0.0;
        q(i,10) = pi/gas_gamma1 + 0.5*(pow(q(i,7),2)+pow(q(i,8),2)
                                      +pow(q(i,9),2))/q(i,6);
        q(i,11) = ex;
        q(i,12) = 0.0;
        q(i,13) = 0.0;
        q(i,14) = 0.0;
        q(i,15) = by;
        q(i,16) = 0.0;
    }
}
