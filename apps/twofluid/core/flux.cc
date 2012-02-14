#include "miniwarpx.h"
#include "twofluid.h"

// flux function for Twofluid equations
void 
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    double rhoe,ue,ve,we,Ee,pe;
    double rhoi,ui,vi,wi,Ei,pi;
    double ex,ey,ez,bx,by,bz;

    Twofluid_Vars *tfv = (Twofluid_Vars*) rd.mvar;
    // gas constant
    double gas_gamma = tfv->gas_gamma;
    double gas_gamma1 = gas_gamma-1;
    // speed of light
    double c0 = tfv->c0;

    for (int i=1-mbc; i<=mx+mbc; i++)
    {
        // compute primitive variables
        
        // electon fluid
        rhoe = q(i,1);
        ue   = q(i,2)/rhoe;
        ve   = q(i,3)/rhoe;
        we   = q(i,4)/rhoe;
        Ee   = q(i,5);
        pe   = gas_gamma1*(Ee-0.5*rhoe*(ue*ue+ve*ve+we*we));
        // ion fluid
        rhoi = q(i,6);
        ui   = q(i,7)/rhoi;
        vi   = q(i,8)/rhoi;
        wi   = q(i,9)/rhoi;
        Ei   = q(i,10);
        pi   = gas_gamma1*(Ei-0.5*rhoi*(ui*ui+vi*vi+wi*wi));
        // electromagnetic field
        ex   = q(i,11);
        ey   = q(i,12);
        ez   = q(i,13);
        bx   = q(i,14);
        by   = q(i,15);
        bz   = q(i,16);
        
        // compute flux

        // electron fluid
        fx(i,1) = rhoe*ue;
        fx(i,2) = rhoe*ue*ue + pe;
        fx(i,3) = rhoe*ue*ve;
        fx(i,4) = rhoe*ue*we;
        fx(i,5) = (Ee+pe)*ue;

        // ion fluid
        fx(i,6) = rhoi*ui;
        fx(i,7) = rhoi*ui*ui + pi;
        fx(i,8) = rhoi*ui*vi;
        fx(i,9) = rhoi*ui*wi;
        fx(i,10) = (Ei+pi)*ui;

        // electromagentic fields
        fx(i,11) = 0;
        fx(i,12) = c0*c0*bz;
        fx(i,13) = -c0*c0*by;
        fx(i,14) = 0;
        fx(i,15) = -ez;
        fx(i,16) = ey;
    }
}
