#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "miniwarpx.h"
#include "farray.h"
#include "euler.h"

/**
   Writes conserved variables data to file. Only the cell averages are
   written out by this version.

   Parameters
   ----------

   rd [in]    - Data for this simulations
   fname [in] - Output file name
   q [in]     - Conserved variables

   Notes
   -----

   The output is an ASCII file. There are 'rd.mx' rows, and in each
   row the coefficents of the basis function are written. This file
   can be read into matlab using a load('frame.xx'); command. The read
   data should be reshaped before processing.
 */
void 
out(const Run_Data& rd, char const *fname, const FArray<double>& q)
{
    Euler_Vars *ev = (Euler_Vars*) rd.mvar;

    char buff[256];
    sprintf(buff, "./%s/%s", rd.run_name, fname);
    // open file to write solution
    std::ofstream fout(buff);

    // loop over all interior points
    for(int i=1; i<=rd.mx; i++)
    {
        for(int m=1; m<=rd.meqn; m++)
            for(int c=1; c<=rd.ncoeffs; c++)
                fout << q(i,m,c) << " "; // write coefficient to file

        fout << ev->ef(i,1) << " "; // write electric field to file

        fout << std::endl;
    }
}
