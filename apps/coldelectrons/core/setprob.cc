#include <iostream>
#include "miniwarpx.h"
#include "cold.h"

using namespace std;

void
setprob(Run_Data& rd)
{

    char *val;
    Cold_Vars *cv = new Cold_Vars();

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

    // c0
    val = lookup_in_section(s, "c0");
    if(val)
        cv->c0 = extract_double(val);
    else
        // assume c0 = 1.0
        cv->c0 = 1.0;

    // electron charge
    val = lookup_in_section(s, "q");
    if(val)
        cv->q = extract_double(val);
    else
        // assume ion charge is 1
        cv->q = 1.0;

    // electron mass
    val = lookup_in_section(s, "m");
    if(val)
        cv->m = extract_double(val);
    else
        // assume ion mass is 1
        cv->m = 1.0;

    // epsilon0
    val = lookup_in_section(s, "epsilon0");
    if(val)
        cv->epsilon0 = extract_double(val);
    else
        // assume epsilon0 is 1
        cv->epsilon0 = 1.0;

    rd.mvar = (void*) cv; // cast is needed to store stuff properly
}
