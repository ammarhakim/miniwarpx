#include "miniwarpx.h"
#include <iostream>

// dummy function: does nothing. The user MUST override this if the
// simulation is to work properly

/*
  Sets the initial condition.

  Parameters
  ----------

  rd [in]   - Input data for simulation
  xloc [in] - xloc(i) is the coordinate at which the conserved variable should be computed
  q [in]    - q(i,m) is the mth conserved variable in cell i computed at xloc(i)

*/
void 
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    std::cout << "** MINIWARPX: qinit is not provided or is not linked properly" << std::endl;
}
