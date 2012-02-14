#include "miniwarpx.h"

/**
   Sets up 'rd' with MACCOR2 algorithm specific data. This essentially
   means allocating memory for the 'maccor2_step' routine.
*/
void
maccor2_setup(Run_Data &rd)
{

    int mx = rd.mx;
    int mbc = rd.mbc;
    int meqn = rd.meqn;

    rd.ncoeffs = 1; // MACCOR2 needs one level of storage
    rd.step = maccor2_step; // function to advance solution by dt
    
    // allocate workspace for MACCOR2 algorithm: the assignment
    // operator rellocates memory for the needed arrays.
    
    MACCOR2_Workspace *mws = new MACCOR2_Workspace();
    mws->s    = FArray<double>(Range(1-mbc,mx+mbc),0.0);
    mws->q    = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    mws->f    = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    mws->fl   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    mws->fr   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    mws->sr   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn),0.0);
    mws->qp   = FArray<double>(Range(1-mbc,mx+mbc, 1,meqn, 1,1),0.0);

    mws->dtdx   = FArray<double>(Range(1-mbc,mx+mbc),0.0);

    rd.work = (void*) mws;
}
