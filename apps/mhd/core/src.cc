#include "mhd.h"
#include "miniwarpx.h"

inline
static
double
pressure(double gas_gamma, const FArray<double>& q, int i)
{
    // rho*(u^2+v^2+w^2)
    double rhou2 = (pow(q(i,2),2) + pow(q(i,3),2) + pow(q(i,4),2))/q(i,1);
    // Bx^2+By^2+Bz^2
    double B2 = pow(q(i,6),2) + pow(q(i,7),2) + pow(q(i,8),2);

    // compute pressure
    return (gas_gamma-1)*(q(i,5) - 0.5*rhou2 - 0.5*B2);
}

void 
src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q)
{
    int mx = rd.mx;

    MHD_Vars *mv = (MHD_Vars*) rd.mvar;
    double gas_gamma = mv->gas_gamma; // gas constant
    double r,rho,u,v,w,pr,bx,by,bz;

    bool is_radial = mv->is_radial;
    if(is_radial==1)
    {
        // loop over all interior cells computing geometric source terms
        for(int i=1; i<=mx; i++)
        {
            r = rd.xcoords(i); // cell center coordinate
            
            // compute primitive variables
            rho = q(i,1);
            u = q(i,2)/rho;
            v = q(i,3)/rho;
            w = q(i,4)/rho;
            
            pr = pressure(gas_gamma,q,i); // fluid pressure
            
            bx = q(i,6);
            by = q(i,7);
            bz = q(i,8);
            
            pr = pr + 0.5*(bx*bx + by*by + bz*bz); // total pressure
            
            // compute source term
            sr(i,1) = -rho*u/r;
            sr(i,2) = -(rho*u*u + by*by)/r;
            sr(i,3) = -2.0*rho*u*v/r;
            sr(i,4) = -rho*u*w/r;
            sr(i,5) = -u*(q(i,5)+pr)/r;
            sr(i,6) = 0.0;
            sr(i,7) = -u*by/r;
            sr(i,8) = 0.0;//(w*bx-u*bz)/r;
        }
    }
}
