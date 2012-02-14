#include "miniwarpx.h"

/**
   Computes fluctuations
 */
void 
eval_fluctuations(Run_Data& rd, const FArray<double>& wave, const FArray<double>& s,
                  FArray<double>& amdq, FArray<double>& apdq)
{
    if(rd.edge_splitting==f_wave)
        eval_fluctuations_fwave(rd,wave,s,amdq,apdq);
    else
        eval_fluctuations_qwave(rd,wave,s,amdq,apdq);
}
