#include <iostream>
#include "miniwarpx.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    // Problem specific initializations go in this function. For the
    // advextion equation the advextion velocity is read and set in
    // 'rd.rpar'.

    char *val;

    // read advection velocity and set it
    rd.nrpar = 1; // one extra double is needed
    rd.rpar = new double[rd.nrpar]; // allocate memory for it

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
    val = lookup_in_section(s, "adv_speed");
    if(val)
        // store advection velocity in 'rpar'
        rd.rpar[0] = extract_double(val);
    else
        // assume velocity is 1
        rd.rpar[0] = 1.0;
}
