#include "miniwarpx.h"

/**
   Computes fluctuations for f-wave method
 */
void 
eval_fluctuations_fwave(Run_Data& rd, const FArray<double>& wave, const FArray<double>& s,
                        FArray<double>& amdq, FArray<double>& apdq)
{
    int mx = rd.mx;
    int mbc = rd.mbc;    
    int meqn = rd.meqn;
    int mwave = rd.mwave;

    for(int i=2-mbc; i<=mx+mbc; i++)
        for(int m=1; m<=meqn; m++)
        {
            amdq(i,m) = 0.0;
            apdq(i,m) = 0.0;
            for(int mw=1; mw<=mwave; mw++)
            {
                if(s(i,mw)<0)
                    // left going wave
                    amdq(i,m) += wave(i,m,mw);
                else if(s(i,mw)>0)
                    // right going wave
                    apdq(i,m) += wave(i,m,mw);
                else
                { // zero wave: add contribution to both apdq and amdq
                    apdq(i,m) += 0.5*wave(i,m,mw);
                    amdq(i,m) += 0.5*wave(i,m,mw);
                }
            }
        }
}
