#include "miniwarpx.h"
#include <iostream>

// functions named qinit and rp must be provided by the user and
// linked properly to create the xminiwarpx executable
extern void qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q);
extern void setkappa(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q);
extern void flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q);
extern void src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q);
extern void maxs(const Run_Data& rd, FArray<double>& s, FArray<double>& q);
extern void rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& f,
               FArray<double>& wave, FArray<double>& s, 
               FArray<double>& amdq, FArray<double>& apdq);

using namespace std;

extern void copy_file(char *s, char *d);

/**
  Reads data from file 'fname' and if inputs are consistent, runs the
  simulation

  Parameters
  ----------

  rd [in/out] - Uninitialized Run_Data object

 */
void
driver(Run_Data &rd, char *fname)
{
    int passed;

    //
    // Parse input file and set simlulation information
    //
    passed = read_inp(rd, fname);
    if(passed == 0)
    { // input file was not parsed properly: abort
        cout << "** MINIWARPX: Unable to parse input file. Aborting..." << endl;
        exit(0);
    }

    rd.inp_file = strdup(fname); // set input file name for use later

    // at this point data consistency should be verified

    // set pointers to various function: users could potentially
    // change these in setprob if so desired
    rd.rp = rp; // Reimann solver
    rd.qinit = qinit; // initial conditions
    rd.flux = flux; // flux function
    rd.src = src; // source function
    rd.maxs = maxs; // maximum wave speed
    rd.setkappa = setkappa; // capacity function

    //
    // Call user specified routine to set application specific data
    //
    setprob(rd);

    // setup directory for writing output
    passed = init_output(rd);
    if (passed != 1)
    { // there was some problem with the directory creation process
        cout << "** MINIWARPX: Unable to initialize output. Aborting... " << endl;
        exit(0);
    }

    // copy the input file to the output directory so that the
    // simulation can be done again later if so desired,
    char buff[256];
    sprintf(buff, "./%s/%s", rd.run_name, rd.inp_file);
    copy_file(fname, buff);

    // 
    // Run the simulation
    //
    solve(rd);
}
