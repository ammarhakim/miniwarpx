#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "maxwell.h"
#include "miniwarpx.h"
#include "farray.h"

/**
   Write epsilon and mu to file
 */
void
write_coeffs(const Run_Data& rd, char const *fname)
{
    IMaxwell_Vars *v = (IMaxwell_Vars*) rd.mvar;
    char buff[256];
    sprintf(buff, "./%s/%s", rd.run_name, fname);
    // open file to write solution
    std::ofstream fout(buff);

    for(int i=1; i<=rd.mx; ++i)
        fout << v->ep(i) << " " << v->mu(i) << std::endl;
}
