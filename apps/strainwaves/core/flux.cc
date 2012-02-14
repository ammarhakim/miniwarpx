#include "miniwarpx.h"
#include "strainwave.h"

//flux function for Maxwell's equations
void
flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    

    Strainwave_Vars *sv = (Strainwave_Vars*) rd.mvar;

    for (int i=2-mbc; i<=mx+mbc; i++)
    {
        // compute flux
        fx(i,1) = -q(i,2)/sv->rho(i);;
        fx(i,2) = -1*stress(sv->beta,sv->modul(i),q(i,1));
    }
}
