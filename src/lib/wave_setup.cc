#include "miniwarpx.h"

/**
   Sets up 'rd' with WAVE algorithm specific data. This essentially
   means allocating memory for the 'wave_step' routine.
*/
void
wave_setup(Run_Data &rd)
{

    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;
    int mwave = rd.mwave;

    rd.ncoeffs = 1; // WAVE needs one level of storage
    rd.step = wave_step; // function to advance solution by dt
    
    // allocate workspace for WAVE algorithm: the assignment
    // operator rellocates memory for the needed arrays.
    
    WAVE_Workspace *wws = new WAVE_Workspace();
    wws->s    = FArray<double>(Range(1-mbc,mx+mbc, 1,mwave),0.0);
    wws->wave = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,mwave),0.0);
    wws->amdq = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    wws->apdq = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    wws->fs   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    wws->f    = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    wws->q    = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    wws->fl   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    wws->fr   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);

    wws->dtdx   = FArray<double>(Range(1-mbc,mx+mbc),0.0);

    rd.work = (void*) wws;

}
