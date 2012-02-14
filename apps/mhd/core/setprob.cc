#include <iostream>
#include "miniwarpx.h"
#include "mhd.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    // Problem specific initializations go in this function. For the
    // MDH equation the gas adiabatic index must be specified

    char *val;
    MHD_Vars *mv = new MHD_Vars();  // allocate memory for MHD_Vars object

    // get a hold of the [APPLICATION-DATA] section
    Section *s;
    s = lookup_in_chapter(rd.inpdata, "APPLICATION-DATA");
    if (s==NULL)
    {
        cout << "Missing [APPLICATION-DATA] in input file " 
             << rd.inp_file
             << endl;
        exit(1); // rather ungracious exit!
    }
    val = lookup_in_section(s, "gas_gamma");
    if(val)
        mv->gas_gamma = extract_double(val);
    else
        // assume it is 1.4
        mv->gas_gamma = 1.4;

    val = lookup_in_section(s, "is_radial");
    if(val)
        mv->is_radial = extract_int(val);
    else
        // assume cartesian
        mv->is_radial = 0;

    rd.mvar = (void*) mv; // cast is needed to store stuff properly
}
