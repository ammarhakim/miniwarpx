#include <iostream>
#include "miniwarpx.h"
#include "twofluid.h"

using namespace std;

void
setprob(Run_Data& rd)
{
    // Problem specific initializations go in this function. For the
    // twofluid equations the gas adiabatic index, speed of light,
    // ion charge and mass, electron charge and mass and epsilon0
    // are read and set in 'rd.mvar'

    char *val;
    Twofluid_Vars *tfv = new Twofluid_Vars(); // allocate memory for Twofluid_Var object

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
        tfv->gas_gamma = extract_double(val);
    else
        // assume gas_gamma is 1.4
        tfv->gas_gamma = 1.4;

    // speed of light, c0
    val = lookup_in_section(s, "c0");
    if(val)
        tfv->c0 = extract_double(val);
    else
        // assume c0 is 1.0
        tfv->c0 = 1.0;

    // ion charge
    val = lookup_in_section(s, "qi");
    if(val)
        tfv->qi = extract_double(val);
    else
        // assume ion charge is 1
        tfv->qi = 1.0;

    // electron charge
    val = lookup_in_section(s, "qe");
    if(val)
        tfv->qe = extract_double(val);
    else
        // assume electron charge is -1
        tfv->qe = -1.0;

    // ion mass
    val = lookup_in_section(s, "mi");
    if(val)
        tfv->mi = extract_double(val);
    else
        // assume ion mass is 1
        tfv->mi = 1.0;

    // electron mass
    val = lookup_in_section(s, "me");
    if(val)
        tfv->me = extract_double(val);
    else
        // assume electron mass is 1
        tfv->me = 1.0;

    // epsilon0
    val = lookup_in_section(s, "epsilon0");
    if(val)
        tfv->epsilon0 = extract_double(val);
    else
        // assume epsilon0 is 1
        tfv->epsilon0 = 1.0;

    val = lookup_in_section(s, "is_radial");
    if(val)
        tfv->is_radial = extract_int(val);
    else
        // assume cartesian
        tfv->is_radial = 0;

    int maxmx = rd.mx + rd.mbc; // for allocating memory

    // allocate memory for work arrays needed in Reimann solver
    tfv->u2v2w2 = FArray<double>(Range(-1,maxmx), 0.);
    tfv->u      = FArray<double>(Range(-1,maxmx), 0.);
    tfv->v      = FArray<double>(Range(-1,maxmx), 0.);
    tfv->w      = FArray<double>(Range(-1,maxmx), 0.);
    tfv->enth   = FArray<double>(Range(-1,maxmx), 0.);
    tfv->a      = FArray<double>(Range(-1,maxmx), 0.);
    tfv->g1a2   = FArray<double>(Range(-1,maxmx), 0.);    
    tfv->euv    = FArray<double>(Range(-1,maxmx), 0.);    

    rd.mvar = (void*) tfv; // cast is needed to store stuff properly
}
