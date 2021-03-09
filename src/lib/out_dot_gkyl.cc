#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>

#include "miniwarpx.h"
#include "farray.h"

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
out_dot_gkyl(const Run_Data& rd, char const *fname, const FArray<double>& q)
{

    char buff[256];
    sprintf(buff, "./%s/%s", rd.run_name, fname);
    // open file to write solution
    FILE *fp = fopen(buff, "wb"); if (!fp) return;

    // write grid
    uint64_t ndim = 1;
    uint64_t cells[] = { rd.mx };
    double lower[] = { rd.xlower }, upper[] = { rd.xupper };

    uint64_t real_type = 2;
    fwrite(&real_type, sizeof(uint64_t), 1, fp);

    fwrite(&ndim, sizeof(uint64_t), 1, fp);
    fwrite(cells, sizeof(uint64_t), 1, fp);
    fwrite(lower, sizeof(double), 1, fp);
    fwrite(upper, sizeof(double), 1, fp);

    uint64_t esznc = sizeof(double)*rd.meqn*rd.ncoeffs, size = rd.mx;
    fwrite(&esznc, sizeof(uint64_t), 1, fp);
    fwrite(&size, sizeof(uint64_t), 1, fp);

    std::vector<double> data(rd.mx*rd.meqn*rd.ncoeffs);
    int count = 0;
    for(int i=1; i<=rd.mx; i++)
    {
        for(int m=1; m<=rd.meqn; m++)
            for(int c=1; c<=rd.ncoeffs; c++)
              data[count++] = q(i,m,c);
    }

    fwrite(&data[0], esznc*size, 1, fp);
    fclose(fp);
}
