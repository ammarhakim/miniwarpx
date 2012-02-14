#include "tenmom.h"
#include "miniwarpx.h"
#include <math.h>

void 
src_p(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{

    double r[3],u[3],v[3],w[3];
    double pxx[3],pxy[3],pxz[3],pyy[3],pyz[3],pzz[3];
    double ex,ey,ez,bx,by,bz;
    double re,ri,epsilon0;

    // set constants
    Tenmom_Vars *tmv = (Tenmom_Vars*) rd.mvar;

    re = tmv->qe/tmv->me; // electron charge-mass ratio
    ri = tmv->qi/tmv->mi; // ion charge-mass ratio
    epsilon0 = tmv->epsilon0; // permittivity of free space

//     compute primitive variables
//    -- electron quantities
    r[1]   = q(1);
    u[1]   = q(2)/q(1);
    v[1]   = q(3)/q(1);
    w[1]   = q(4)/q(1);
    pxx[1] = q(5);
    pxy[1] = q(6);
    pxz[1] = q(7);
    pyy[1] = q(8);
    pyz[1] = q(9);
    pzz[1] = q(10);
//     -- ion quantities
    r[2]   = q(11);
    u[2]   = q(12)/q(11);
    v[2]   = q(13)/q(11);
    w[2]   = q(14)/q(11);
    pxx[2] = q(15);
    pxy[2] = q(16);
    pxz[2] = q(17);
    pyy[2] = q(18);
    pyz[2] = q(19);
    pzz[2] = q(20);
//     -- electromagnetic quantities
    ex = q(21);
    ey = q(22);
    ez = q(23);
    bx = q(24);
    by = q(25);
    bz = q(26);

//     compute electron source terms
//     -- mass source term
    sr(1) = 0.0;
//     -- momentum source terms
    sr(2) = r[1]*re*(ex + v[1]*bz - w[1]*by);
    sr(3) = r[1]*re*(ey + w[1]*bx - u[1]*bz);
    sr(4) = r[1]*re*(ez + u[1]*by - v[1]*bx);
//     -- pressure tensor source terms
    sr(5) = 2.0*r[1]*re*u[1]*ex  +  2.0*re*(bz*pxy[1] - by*pxz[1]);
    sr(6) = r[1]*re*(u[1]*ey + v[1]*ex) +
        re*(bz*pyy[1] - by*pyz[1] - bz*pxx[1] + bx*pxz[1]);
    sr(7) = r[1]*re*(u[1]*ez + w[1]*ex)  + 
        re*(bz*pyz[1] + by*pxx[1] - by*pzz[1] - bx*pxy[1]);
    sr(8) = 2.0*r[1]*re*v[1]*ey  +  2.0*re*(bx*pyz[1] - bz*pxy[1]);
    sr(9) = r[1]*re*(v[1]*ez + w[1]*ey)  + 
        re*(by*pxy[1] - bz*pxz[1] + bx*pzz[1] - bx*pyy[1]);
    sr(10) = 2.0*r[1]*re*w[1]*ez  +  2.0*re*(by*pxz[1] - bx*pyz[1]);

//     compute ion source terms
//     -- mass source term
    sr(11) = 0.0;
//     -- momentum source terms
    sr(12) = r[2]*ri*(ex + v[2]*bz - w[2]*by);
    sr(13) = r[2]*ri*(ey + w[2]*bx - u[2]*bz);
    sr(14) = r[2]*ri*(ez + u[2]*by - v[2]*bx);
//     -- pressure tensor source terms
    sr(15) = 2.0*r[2]*ri*u[2]*ex  +  2.0*ri*(bz*pxy[2] - by*pxz[2]);
    sr(16) = r[2]*ri*(u[2]*ey + v[2]*ex) +
        ri*(bz*pyy[2] - by*pyz[2] - bz*pxx[2] + bx*pxz[2]);
    sr(17) = r[2]*ri*(u[2]*ez + w[2]*ex) +
        ri*(bz*pyz[2] + by*pxx[2] - by*pzz[2] - bx*pxy[2]);
    sr(18) = 2.0*r[2]*ri*v[2]*ey  +  2.0*ri*(bx*pyz[2] - bz*pxy[2]);
    sr(19) = r[2]*ri*(v[2]*ez + w[2]*ey)  + 
        ri*(by*pxy[2] - bz*pxz[2] + bx*pzz[2] - bx*pyy[2]);
    sr(20) = 2.0*r[2]*ri*w[2]*ez  +  2.0*ri*(by*pxz[2] - bx*pyz[2]);

//     compute electromagnetic source terms
//    -- electric field source terms
    sr(21) = -(re*r[1]*u[1] + ri*r[2]*u[2])/epsilon0;
    sr(22) = -(re*r[1]*v[1] + ri*r[2]*v[2])/epsilon0;
    sr(23) = -(re*r[1]*w[1] + ri*r[2]*w[2])/epsilon0;
//     -- magnetic field source terms
    sr(24) = 0.0;
    sr(25) = 0.0;
    sr(26) = 0.0;
}
