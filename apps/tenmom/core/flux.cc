#include <math.h>
#include "miniwarpx.h"

// flux function for Tenmoment equations
void
flux(const Run_Data& rd, FArray<double>& f, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;
    
    double r,u,v,w,pxx,pxy,pxz,pyy,pyz,pzz;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        r = q(i,1);
        u = q(i,2)/q(i,1);
        v = q(i,3)/q(i,1);
        w = q(i,4)/q(i,1);
        pxx = q(i,5)-r*u*u;
        pxy = q(i,6)-r*u*v;
        pxz = q(i,7)-r*u*w;
        pyy = q(i,8)-r*v*v;
        pyz = q(i,9)-r*v*w;
        pzz = q(i,10)-r*w*w;
        
        // compute flux
        f(i,1) = r*u;
        f(i,2) = r*u*u + pxx;
        f(i,3) = r*u*v + pxy;
        f(i,4) = r*u*w + pxz;
        f(i,5) = r*u*u*u + 3.0*u*pxx;
        f(i,6) = r*u*u*v + 2.0*u*pxy + v*pxx;
        f(i,7) = r*u*u*w + 2.0*u*pxz + w*pxx;
        f(i,8) = r*u*v*v + u*pyy + 2.0*v*pxy;
        f(i,9) = r*u*v*w + u*pyz + v*pxz + w*pxy;
        f(i,10) = r*u*w*w + u*pzz + 2.0*w*pxz;

    }
}
