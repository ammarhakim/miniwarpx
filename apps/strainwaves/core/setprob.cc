#include <iostream>
#include "miniwarpx.h"
#include "strainwave.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    // Problem specific initializations go in this function. For the
    // Maxwell equation the speed of light is read and set in
    // 'rd.rpar'.

    char *val;
    Strainwave_Vars *sv = new Strainwave_Vars();  // allocate memory for Strainwave_Vars object

    sv->modul  = FArray<double>(Range(1-rd.mbc,rd.mx+rd.mbc), 0.);
    sv->rho    = FArray<double>(Range(1-rd.mbc,rd.mx+rd.mbc), 0.);

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
    val = lookup_in_section(s, "ubar");
    if(val)
        // store initial speed in 'rpar'
        sv->ubar = extract_double(val);
    else
        // assume initial speed is 1
        sv->ubar = 1.0;

    val = lookup_in_section(s, "beta");
    if(val)
        // store linear/nonlinear in 'rpar'
        sv->beta = extract_double(val);
    else
        // assume linear
        sv->beta = 0.0;

    rd.mvar = (void*) sv; // cast is needed to store stuff properly
}
