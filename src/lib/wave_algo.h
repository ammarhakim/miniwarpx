#ifndef __wave_algo__2005__
#define __wave_algo__2005__

#include "miniwarpx.h"

//
// This file is part of WarpX hyperbolic conservation law library.
//
// Author: Ammar Hakim
// Date  : November 2005
//

/*
  Data and functions specific to the WAVE algorithm
 */

void wave_source_advance(Run_Data& rd, FArray<double>& q,double tcurr, double dt);

#endif // __wave_algo__2005__
