#include <stdio.h>
#include <math.h>
#include "miniwarpx.h"
#include "strainwave.h"

void
qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q)
{
    int num = 0;
    double xcell;
    double pi   = 3.14159;
    double del  = 1.0;
    double alpha= 0.5;        // for layering

    // initial conditions
    double rhoa = 4.0;          // density of material 'a'
    double rhob = 1.0;          // density of material 'b'
    double Ka   = 4.0;       // bulk modulus of material 'a'
    double Kb   = 1.0;          // bulk modulus of material 'b'

    Strainwave_Vars *sv = (Strainwave_Vars*) rd.mvar;


     for(int i=2-rd.mbc; i<=rd.mx+rd.mbc; i++)
     {
//         for fewer layers use this part
//         xcell = xloc(i);
//         if (xcell<=15 || xcell>=35)
//         {
//             sv->rho(i) = rhoa;
//             sv->modul(i) = Ka;
//          }
//         else if (xcell>15 && xcell<35)
//          {
//              sv->rho(i)   = rhob;
//              sv->modul(i) = Kb;
//          }

//     this part is for setting the layered medium so that at 
//     every 0.5 in x, the medium alternates

        xcell = xloc(i);

        if ( num*del<=xcell && xcell<(num+alpha)*del)
        {
            sv->rho(i)   = rhoa;
            sv->modul(i) = Ka;
        }
        else
        {
            sv->rho(i)   = rhob;
            sv->modul(i) = Kb;
        }
        if (xcell>=num+1-0.05) num++;

        q(i,1) = 0.0;
        q(i,2) = 0.0;
        
    }
}
