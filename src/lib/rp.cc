#include "miniwarpx.h"
#include <iostream>

// dummy function: does nothing. The user MUST override this if the
// simulation is to work properly


/*
  Reimann solver for equation system.

  Parameters
  ----------
  
  rd [in] - Input data for simulation
  ql [in] - Conserved variables evaluated at left cell edge. Thus ql(i,*) is the 
            conserved variable at the left edge of cell i
  qr [in] - Conserved variables evaluated at right cell edge. Thus qr(i,*) is the 
            conserved variable at the right edge of cell i
  df [in] - Vector to split at each cell edge. This is usuall the jump in q at each
            edge or the jump in flux at each edge. However, for the RKDG algorithm it
            could be some other vector too. The user should not make any assumptions
            about it in this routine
  wave [out] - wave(i,m,mw) is the mth component of wave mw at edge i.
  s    [out] - s(i,mw) is the speed of wave mw
  amdq [out] - amdq(i,m) is the mth component of the left-going fluctuation at edge i
  apdq [out] - apdq(i,m) is the mth component of the right-going fluctuation at edge i


  Notes
  -----

  For definition of waves and fluctuations see LeVeque's book or WarpX
  notes. For a cell i its left edge is labeled i and its right edge is
  labelled i+1. Thus there is one extra cell edge than cells. 

  To compute the fluctuations, in most cases the provided routine
  eval_fluctuations() can be used. This automatically takes care of
  the slight difference in the definition of the fluctuations
  depending on if the jump in q (q_waves) or the jump in flux (f_wave)
  approach is being used.

*/
void 
rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
   FArray<double>& wave, FArray<double>& s, 
   FArray<double>& amdq, FArray<double>& apdq)
{
    std::cout << "** MINIWARPX: rp is not provided or is not linked properly" << std::endl;
}
