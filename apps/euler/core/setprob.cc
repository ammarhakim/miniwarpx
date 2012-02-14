#include <iostream>
#include "miniwarpx.h"
#include "euler.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    // Problem specific initializations go in this function. For the
    // Euler equation the gas adiabatic index is read and set in
    // 'rd.mvar'

    char *val;
    Euler_Vars *ev = new Euler_Vars(); // allocate memory for Euler_Var object

    // read gas_gamma and set it

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
    // gas_gamma
    val = lookup_in_section(s, "gas_gamma");
    if(val)
        ev->gas_gamma = extract_double(val);
    else
        // assume gas_gamma is 1.4
        ev->gas_gamma = 1.4;

    int maxmx = rd.mx + rd.mbc; // for allocating memory

    // allocate memory for work arrays needed in Reimann solver
    ev->u2v2w2 = FArray<double>(Range(-1,maxmx), 0.);
    ev->u      = FArray<double>(Range(-1,maxmx), 0.);
    ev->v      = FArray<double>(Range(-1,maxmx), 0.);
    ev->w      = FArray<double>(Range(-1,maxmx), 0.);
    ev->enth   = FArray<double>(Range(-1,maxmx), 0.);
    ev->a      = FArray<double>(Range(-1,maxmx), 0.);
    ev->g1a2   = FArray<double>(Range(-1,maxmx), 0.);    
    ev->euv    = FArray<double>(Range(-1,maxmx), 0.);    

    rd.mvar = (void*) ev; // cast is needed to store stuff properly
}
