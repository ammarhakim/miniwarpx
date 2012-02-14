#include "miniwarpx.h"
#include <iostream>
#include <string.h>

using namespace std;

/*
  Reads data from the input file 'fname', populating simulation data
  in 'rd'. A pointer to the read data is also stored as 'rd.inpdata'
  so that applications can initialize application specific data later
  if desired.

  Parameters
  ----------

  rd    [in/out] - Uninitialized Run_Data object
  fname [in]     - Name of input file

  Returns
  -------

  If the read succeeds returns 1, or 0 otherwise.

 */
int
read_inp(Run_Data& rd, char *fname)
{

    Section *s;
    char *val;
    int passed = 1; // flag to indicate if populating 'rd' worked

    // parse the input data and store it in 'rd'
    rd.inpdata = parse_input(fname);

    if (rd.inpdata==NULL)
    {
        cout << "Parsing input file " << fname << " failed" << endl;
        return 0;
    }

    // Read data from section [RUN-DATA]
    s = lookup_in_chapter(rd.inpdata, "RUN-DATA");
    if (s==NULL)
    {
        cout << "Missing [RUN-DATA] in file " << fname << endl;
        return 0;
    }
    // description
    val = lookup_in_section(s, "description");
    if(val)
        rd.description = strdup(val);
    else
        rd.description = strdup("MINIWARPX: Description not given");

    // run_name
    val = lookup_in_section(s, "run_name");
    if(val)
        rd.run_name = strdup(val);
    else
        rd.run_name = "run";

    // meqn
    val = lookup_in_section(s, "meqn");
    if(val)
        rd.meqn = extract_int(val);
    else
        rd.meqn = 1;

    // mx
    val = lookup_in_section(s, "mx");
    if(val)
        rd.mx = extract_int(val);
    else
        rd.mx = 100;

    // mwave
    val = lookup_in_section(s, "mwave");
    if(val)
        rd.mwave = extract_int(val);
    else
    {
        cout << "No of waves not specified in file " << fname << endl;
        passed = 0;
    }

    // mbc
    val = lookup_in_section(s, "mbc");
    if(val)
        rd.mbc = extract_int(val);
    else
        rd.mbc = 1;

    // nout
    val = lookup_in_section(s, "nout");
    if(val)
        rd.nout = extract_int(val);
    else
        rd.nout = 1;

    // verbose
    val = lookup_in_section(s, "verbose");
    if(val)
        if (strcmp(val, "true")==0)
            rd.verbose = true;
        else
            rd.verbose = false;
    else
        rd.verbose = false;

    // algo
    val = lookup_in_section(s, "algo");
    if(val)
    {
        if (strcmp(val, "WAVE")==0)
            rd.algo = WAVE;
        else if (strcmp(val, "RKDG")==0)
            rd.algo = RKDG;
        else if (strcmp(val, "MACCOR2")==0)
            rd.algo = MACCOR2;
        else
        {
            cout << "Unknown algorithm option " << val << " in file " << fname << endl;
            passed = 0; // the parsing has failed
        }
    }
    else
        rd.algo = WAVE;

    // has_source
    val = lookup_in_section(s, "has_source");
    if(val)
        if (strcmp(val, "true")==0)
            rd.has_source = true;
        else
            rd.has_source = false;
    else
        rd.has_source = false;

    // has_kappa
    val = lookup_in_section(s, "has_kappa");
    if(val)
        if (strcmp(val, "true")==0)
            rd.has_kappa = true;
        else
            rd.has_kappa = false;
    else
        rd.has_kappa = false;

    // edge_splitting
    val = lookup_in_section(s, "edge_splitting");
    if(val)
        if(strcmp(val,"q_wave")==0)
            rd.edge_splitting = q_wave;
        else
            rd.edge_splitting = f_wave;
    else
        rd.edge_splitting = f_wave; // by default use F-wave splitting

    // xlower
    val = lookup_in_section(s, "xlower");
    if(val)
        rd.xlower = extract_double(val);
    else
    {
        cout << "xlower not specified in file " << fname << endl;
        passed = 0;
    }

    // xupper
    val = lookup_in_section(s, "xupper");
    if(val)
        rd.xupper = extract_double(val);
    else
    {
        cout << "xupper not specified in file " << fname << endl;
        passed = 0;
    }

    // tstart
    val = lookup_in_section(s, "tstart");
    if(val)
        rd.tstart = extract_double(val);
    else
    {
        cout << "tstart not specified in file " << fname << endl;
        passed = 0;
    }
    // tend
    val = lookup_in_section(s, "tend");
    if(val)
        rd.tend = extract_double(val);
    else
    {
        cout << "tend not specified in file " << fname << endl;
        passed = 0;
    }
    // dt
    val = lookup_in_section(s, "dt");
    if(val)
        rd.dt = extract_double(val);
    else
        rd.dt = 1.0;

    // cfl
    val = lookup_in_section(s, "cfl");
    if(val)
        rd.cfl = extract_double(val);
    else
    {
        cout << "cfl not specified in file " << fname << endl;
        passed = 0;
    }
    // cflm
    val = lookup_in_section(s, "cflm");
    if(val)
        rd.cflm = extract_double(val);
    else
    {
        cout << "cflm not specified in file " << fname << endl;
        passed = 0;
    }

    // bc_left
    val = lookup_in_section(s, "bc_left");
    if(val)
    {
        if(strcmp(val, "bc_periodic")==0)
            rd.bc_type[0] = bc_periodic;
        else if(strcmp(val, "bc_copy")==0)
            rd.bc_type[0] = bc_copy;
        else if(strcmp(val, "bc_wall")==0)
            rd.bc_type[0] = bc_wall;
        else if(strcmp(val, "bc_custom")==0)
            rd.bc_type[0] = bc_custom;
        else
            rd.bc_type[0] = bc_copy;
    }
    else
    {
        cout << "Left boundary condition not specified in file " << fname << endl;
        passed = 0;
    }
    // bc_right
    val = lookup_in_section(s, "bc_right");
    if(val)
    {
        if(strcmp(val, "bc_periodic")==0)
            rd.bc_type[1] = bc_periodic;
        else if(strcmp(val, "bc_copy")==0)
            rd.bc_type[1] = bc_copy;
        else if(strcmp(val, "bc_wall")==0)
            rd.bc_type[1] = bc_wall;
        else if(strcmp(val, "bc_custom")==0)
            rd.bc_type[1] = bc_custom;
        else
            rd.bc_type[1] = bc_copy;
    }
    else
    {
        cout << "Right boundary condition not specified in file " << fname << endl;
        passed = 0;
    }

    // nipar
    val = lookup_in_section(s, "nipar");
    if(val)
        rd.nipar = extract_int(val);
    else
        rd.nipar = 0;
    // nrpar
    val = lookup_in_section(s, "nrpar");
    if(val)
        rd.nrpar = extract_int(val);
    else
        rd.nrpar = 0;
    // nvars
    val = lookup_in_section(s, "nvars");
    if(val)
        rd.nvars = extract_int(val);
    else
        rd.nvars = 0;
    
    // read data from section [WAVE] or [RKDG]
    if (rd.algo == WAVE)
    { // WAVE algorithm is being used
        s = lookup_in_chapter(rd.inpdata, "WAVE");
        if(s == NULL)
        {
            cout << "Missing [WAVE] in file  " << fname << endl;
            return 0;
        }

        // source_splitting
        val = lookup_in_section(s, "source_splitting");
        if(val)
            rd.source_splitting = extract_int(val);
        else
            rd.source_splitting = 0;

        // wv_order
        val = lookup_in_section(s, "wv_order");
        if(val)
            rd.wv_order = extract_int(val);
        else
            rd.wv_order = 2; // by default use Lax-Wendroff method
        
        // limiters
        val = lookup_in_section(s, "limiters");
        if(val)
            rd.limiters = extract_int_p(val);
        else
        {
            cout << "Limiters not specified in file " << fname << endl;
            passed = 0;
        }
        // set ncoeffs properly
        rd.ncoeffs = 1;
    }
    else if(rd.algo == RKDG)
    { // RKDG algorithm is being used
        s = lookup_in_chapter(rd.inpdata, "RKDG");
        if(s == NULL)
        {
            cout << "Missing [RKDG] in file  " << fname << endl;
            return 0;
        }

        // rk_order
        val = lookup_in_section(s, "rk_order");
        if(val)
            rd.rk_order = extract_int(val);
        else
        {
            cout << "Runge-Kutta time stepping order not specified in file " << fname 
                 << endl;
            passed = 0;
        }

        // sp_order
        val = lookup_in_section(s, "sp_order");
        if(val)
        {
            rd.sp_order = extract_int(val);
            // also set ncoeffs properly
            rd.ncoeffs = rd.sp_order;
        }
        else
        {
            cout << "Spatial descritization order not specified in file " << fname
                 << endl;
            passed = 0;
        }

        // limiters
        val = lookup_in_section(s, "dg_limiters");
        if(val)
            rd.dg_limiters = extract_int(val);
        else
            rd.dg_limiters = 0; // by default do not apply any limiters to RKDG

        // mmM  for use in the modified min-mod limiter
        val = lookup_in_section(s, "mmM");
        if(val)
            rd.mmM = extract_double(val);
        else
            rd.mmM = 0.0; // by default set it to 0.0
    }   
    else if(rd.algo == MACCOR2)
    { // MacCormick algorithm is being used: no special variables need to be set
    }

    return passed;
}
