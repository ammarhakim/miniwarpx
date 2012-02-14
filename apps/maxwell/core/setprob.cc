#include <iostream>
#include "miniwarpx.h"
#include "maxwell.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    // Problem specific initializations go in this function. For the
    // Maxwell equation the speed of light is read and set in
    // 'rd.rpar'.

    char *val;
    Maxwell_Vars *mv = new Maxwell_Vars();  // allocate memory for Maxwell_Vars object

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
    val = lookup_in_section(s, "c0");
    if(val)
        // store speed of light in 'rpar'
        mv->c0 = extract_double(val);
    else
        // assume speed of light is 1
        mv->c0 = 1.0;

    rd.mvar = (void*) mv; // cast is needed to store stuff properly
}
