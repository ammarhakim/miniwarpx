#include <iostream>
#include "miniwarpx.h"
#include "maxwell.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    IMaxwell_Vars *v = new IMaxwell_Vars();

    v->ep = FArray<double>(Range(1-rd.mbc,rd.mx+rd.mbc), 0.);
    v->mu = FArray<double>(Range(1-rd.mbc,rd.mx+rd.mbc), 0.);

    rd.mvar = (void*) v; // cast is needed to store stuff properly
}
