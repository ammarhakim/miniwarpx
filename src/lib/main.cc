#include "miniwarpx.h"
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <iostream>

/**
   Main entry point into the miniwarpx code. This routine simply
   parses the command line arguments, getting the input file name, and
   calls 'driver' to do the actual work. In the executable is run
   without any arguments the input file name is assumed to be warpx.inp

   Example usage is

   ./xminiwarpx -f myinputs.inp

   Command line flags
   ------------------

   [-f <filename>] - Name of input file for simulation. Default is
                     warpx.inp

 */
int
main(int argc, char **argv)
{
    Run_Data rd;
    char c, *fname;
    clock_t tstart, tend;

    tstart = clock(); // get time at start of simulation

    fname = strdup("warpx.inp"); // default file name is warpx.inp

    // read command line parameters using the GNU getopt library
    opterr = 0;
    while((c = getopt(argc, argv, "hf:")) != -1)
        switch (c)
        {
            case 'f':
                // file name flag
                free(fname); // free old memory
                fname = strdup(optarg);
                break;

            case 'h':
                // print help message and return
                printf("Usage: xminiwarpx [-f <filename>]\n");
                printf("  <filename> is the name of input file. Defaults to warpx.inp\n\n");

                return 1;

            case '?':
                // incorrect option passed
                if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n",
                             optopt);
                return 1;

            default:
                abort ();
        }

    driver(rd, fname);

    tend = clock(); // get time at end of simulation
    // print time to run the simulation
    std::cout << "... simulation finished in " 
              << (double)(tend-tstart)/(double)CLOCKS_PER_SEC
              << " seconds" 
              << std::endl;
}
