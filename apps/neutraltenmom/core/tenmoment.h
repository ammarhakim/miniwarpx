#ifndef __tenmoment_h__
#define __tenmoment_h__

#include "farray.h"

// data-structure to hold variables for the Ten Moment equation set
struct TenMoment_Vars
{
    double qtomi; //charge to mass ratio for ions
    double qtome; //charge to mass ratio for electrons
    double mi;    //mass of the ions
    double me;    //mass of the electrons
    double qi;    //charge of the ions
    double qe;    //charge of the electrons;
    double c0;    //speed of light
    double muo;   //permeability
};

#endif // __tenmoment_h__
