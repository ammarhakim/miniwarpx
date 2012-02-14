#include "miniwarpx.h"

/**
   Computes fluctuations for q-wave method
 */
void 
eval_fluctuations_qwave(Run_Data& rd, const FArray<double>& wave, const FArray<double>& s,
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
                    amdq(i,m) += s(i,mw)*wave(i,m,mw);
                else
                    // right going wave
                    apdq(i,m) += s(i,mw)*wave(i,m,mw);
            }
        }
}
