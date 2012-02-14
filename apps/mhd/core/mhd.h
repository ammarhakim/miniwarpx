#ifndef __mhd_h__
#define __mhd_h__

// data-structure to hold variables for MHD equation application
struct MHD_Vars
{
    double gas_gamma; // Gas adiabatic index
    int is_radial; // 1 if doing RZ plain problem
};

#endif // __mhd_h__
