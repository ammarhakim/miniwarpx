#include "miniwarpx.h"
#include "hallmhd.h"

// flux function for HallMHD equations
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    
    double dx = rd.dx;

    double ne,Ee;
    double jx, gradpex, B;
    double ni,ui,vi,wi,Ei,pi;
    double ex,ey,ez,bx,by,bz;

    Hallmhd_Vars *hv = (Hallmhd_Vars*) rd.mvar;
    // gas constant
    double gas_gamma = hv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    // speed of light
    double c0 = hv->c0;
    double ur = hv->ur;
    double rli= hv->rli;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables

        // ion fluid
        ni   = q(i,1);
        ui   = q(i,2)/ni;
        vi   = q(i,3)/ni;
        wi   = q(i,4)/ni;
        Ei   = q(i,5);
        pi   = gas_gamma1*(Ei-0.5*ni*(ui*ui+vi*vi+wi*wi));
        // electon fluid
        ne   = q(i,6);
        Ee   = q(i,7);
        // electromagnetic field
        ex   = q(i,8);
        ey   = q(i,9);
        ez   = q(i,10);
        bx   = q(i,11);
        by   = q(i,12);
        bz   = q(i,13);

        // compute grad Pex (perp to B) using finite differencing to get 
        // the cell edge values 
        if (i<=mx+mbc-1) gradpex = gas_gamma1*(q(i+1,7)-q(i,7))/dx;
        else gradpex = gas_gamma1*(q(i,7)-q(i-1,7))/dx;

        // compute B^2 and jx (jx in perp direction to B)
        B   = sqrt(bx*bx + by*by + bz*bz);
        if (B != 0.) jx = -(rli*gradpex/B + ne*ex/B - ni*ui);
        else jx = ni*ui;
        
        // compute flux

        // ion fluid
        fx(i,1) = ni*ui;
        fx(i,2) = ni*ui*ui + pi;
        fx(i,3) = ni*ui*vi;
        fx(i,4) = ni*ui*wi;
        fx(i,5) = (Ei+pi)*ui;

        // electron fluid
        fx(i,6) = ni*ui-jx;
        fx(i,7) = (ni*ui/ne-jx/ne)*(1+gas_gamma1)*Ee;

        // electromagentic fields
        fx(i,8) = 0;
        fx(i,9) = c0*c0/(ur*ur)*bz;
        fx(i,10) = -(c0*c0)/(ur*ur)*by;
        fx(i,11) = 0;
        fx(i,12) = -ez;
        fx(i,13) = ey;
    }
}
