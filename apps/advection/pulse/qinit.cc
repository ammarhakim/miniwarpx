#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{

    double xcell;
    for(int i=1; i<=rd.mx; i++)
    {
      double box = 0.0;
      xcell = xloc(i);
      if (xcell>0.7 && xcell<0.9)
        box = 1.0;
      q(i,1) = exp(-100.0*(xcell-0.25)*(xcell-0.25)) + box;
    }
}
