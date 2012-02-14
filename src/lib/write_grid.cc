#include <iostream>
#include <fstream>
#include "miniwarpx.h"

void 
write_grid(const Run_Data& rd, const FArray<double>& xp)
{
    char buff[256];
    sprintf(buff, "./%s/frame.x", rd.run_name);
    // open file to write solution
    std::ofstream fout(buff);

    // loop and write cell-center coordinates for all cells
    for(int i=1; i<=rd.mx; i++)
        fout << xp(i) << std::endl;
}
