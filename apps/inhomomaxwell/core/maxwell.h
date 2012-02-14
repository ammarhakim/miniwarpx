#ifndef __i_maxwell__
#define __i_maxwell__

#include "farray.h"
#include "miniwarpx.h"

struct IMaxwell_Vars
{
    FArray<double> ep; // epsilon(x)
    FArray<double> mu; // mu(x)
};

void write_coeffs(const Run_Data& rd, char const *fname);

#endif //  __maxwell__
