#ifndef __miniwarpx__2005__
#define __miniwarpx__2005__

//
// This file is part of WarpX hyperbolic conservation law library.
//
// Author: Ammar Hakim
// Date  : November 2005
//

#include "farray.h"
#include "utils.h"

#include "inpparse.h"
#include "extractors.h"

// boundary condition flags
const int bc_periodic = 0;
const int bc_copy = 1;
const int bc_wall = 2;
const int bc_custom = 11;

// algorithm type
const int RKDG = 0; // Runge-Kutta Discontinous Galerkim
const int WAVE = 2; // Wave Propagation
const int MACCOR2 = 4; // MacCormick 2nd order space 2nd time 
const int MACCOR4 = 6; // MacCormick 4th order space 2nd time

// type of edge splitting
const int q_wave = 0;
const int f_wave = 1;

class Run_Data; // forward declaration 

/*
  Typedefs for various function pointers used in the system.
*/

//void qinit(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& q);
typedef void (*qptr_t)(const Run_Data&, const FArray<double>& xloc, FArray<double>&);
//void setkappa(const Run_Data& rd, const FArray<double>& xloc, FArray<double>& kappa);
typedef void (*kptr_t)(const Run_Data&, const FArray<double>& xloc, FArray<double>&);
//void flux(const Run_Data& rd, FArray<double>& fx, FArray<double>& q);
typedef void (*fptr_t)(const Run_Data&, FArray<double>&, FArray<double>&);
//void src(const Run_Data& rd, FArray<double>& sr, FArray<double>& q);
typedef void (*sptr_t)(const Run_Data&, FArray<double>&, FArray<double>&);
//void step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla);
typedef void (*stepf_t)(Run_Data&, FArray<double>&, double, double, double&);
//void rp(Run_Data& rd, FArray<double>& ql, FArray<double>& qr, FArray<double>& df,
//        FArray<double>& wave, FArray<double>& s, 
//        FArray<double>& amdq, FArray<double>& apdq)
typedef void (*rpf_t)(Run_Data&, FArray<double>&, FArray<double>&, FArray<double>& df,
                      FArray<double>&, FArray<double>&, 
                      FArray<double>&, FArray<double>&);
//void (*maxs_t)(const Run_Data& rd, FArray<double>& s, FArray<double>& q);
typedef void (*maxs_t)(const Run_Data&, FArray<double>&, FArray<double>&);

// Classes to hold workspace arrays for each method

// Data for WAVE
class WAVE_Workspace
{
public:
    FArray<double> s,wave,amdq,apdq; // speeds,waves,(-/+)fluctuations
    FArray<double> fl,fr; // fluxes adjacent to interfaces
    FArray<double> f; // fluxes
    FArray<double> q; // temporary to hold conserved variables
    FArray<double> fs; // second order correction fluxes
    FArray<double> dtdx; // ratio dt/(dx*kappa)
};

// Data for RKDG
class RKDG_Workspace
{
public:
    FArray<double> w, x; // weights and abscissa for Gaussian quadrature
    FArray<double> pmx, ppmx; // Legendre polynomials and its derivative at x[]
    FArray<double> Cconst; // nomalization constants for basis function
    FArray<double> s,wave,amdq,apdq; // speeds,waves,(-/+)fluctuations
    FArray<double> ql,qr; // conserved variables at right and left edges
    FArray<double> fl,fr; // fluxes adjacent to interfaces
    FArray<double> f;  // fluxes adjacent to interfaces
    FArray<double> df; // jump in flux across an interface
    FArray<double> fedge; // interface fluxes
    FArray<double> q1; // temporary storage for Runge-Kutta time stepping
    FArray<double> rhs; // storage for RHS in Runge-Kutta time stepping

    FArray<double> wave1,wave2; // for use in limiters
};

// Data for MACCOR2
class MACCOR2_Workspace
{
public:
    FArray<double> s; // maxmium wave speed
    FArray<double> f,fl,fr; // fluxes
    FArray<double> sr; // source terms
    FArray<double> q; // temporary storage for conserved variables
    FArray<double> qp; // predicted solution
    FArray<double> dtdx; // ratio dt/(dx*kappa)
};

// Main class specifing simulation inputs
class Run_Data
{
public:

    char *description; // description of this simulation
    char *run_name; // name of run (no spaces allowed)
    // A sub-directory called 'run_name' is created in the directory
    // in which xminiwarpx is executed. All output data from the
    // simulation is written in it.
    char *inp_file; // name of input file

    int meqn; // no of equations in system
    int mx; // no of cells in domain
    int mwave; // no of waves in system
    int mbc; // no ghost cells
    int nout; // no of output files to write
    int verbose; // set to true if more diagonostics are needed

    int algo; // algorithm to use: either RKDG or WAVE

    int has_source; // true if source term is present, false otherwise

    int has_kappa; // true if capacity function is present

    int edge_splitting; // one of q_wave or f_wave
    // for q_waves the jump in the conserved variables across the
    // interface is split while for f_waves the jump in flux across
    // the interface is split

    double xlower; // x coordinate of left edge
    double xupper; // x coordinate of right edge
    double tstart, tend; // start and ending time of simulation
    double cfl; // desired CFL number
    double cflm; // maximum CFL number. cflm>cfl
    double dt; // initial dt to use

    int bc_type[2]; // boundary condition type at left and right edges

    // function pointer for initial conditions
    qptr_t qinit;
    // function pointer for capacity function
    kptr_t setkappa;
    // function pointer for flux function
    fptr_t flux;
    // function pointer for source function
    sptr_t src;
    // function pointer for Reimann solver for the system
    rpf_t rp;
    // function pointer for maximum wave speed calculator
    maxs_t maxs;

    //**
    // RK-DG specific data
    int rk_order; // order of RK Scheme: either 2nd or 3rd order
    int sp_order; // order of spatial scheme: either 2nd or 3rd order
    int dg_limiters; // limiter to apply: one of 
    // 0 - no limiter
    // 1 - characteristics based TVDM limiter
    // 2 - component based non-TVDM limiter: this is much faster but not very good
    double mmM; // value of M for use in minmod function. Usually this is 0.0
    //**

    //**
    // Wave-Propagation specific data
    int source_splitting; // 1 - Gudonov splitting, 2 - Strang splitting.
    int *limiters; // type of limiter for each wave. Should be one of
    // 0 - no limiter applied
    // 1 - minmod
    // 2 - superbee
    // 3 - van Leer
    // 4 - monotonized centered
    // 5 - Beam-Warming
    //**
    int wv_order; // wave splitting method: one of 
    // 1 - Gudonov method
    // 2 - Lax-Wendroff method.
    // Lax-Wendroff is second order accurate and hence to be preferred

    // following two variables communicate values to the user's flux
    // and source routines.
    double tcurrent; // the current time at which flux/src is being called
    FArray<double> xcoords; // coordinate at which flux/src are being evaluated
    // tcurrent and xcoords are useful when the flux or the source
    // depends explicitly on time and/or position.

    // following are provided so that the user can store problem
    // specific variables for use in his supplied routines.

    int nipar; // no of extra integer parameters
    int *ipar; // array of integer parameters: this must be allocated by user

    int nrpar; // no of extra double parameters
    double *rpar; // array of double parameters: this must be allocated by user

    void *mvar; // pointer to user specified data: this must be allocated by user

    int nvars; // no of problem specific variables
    void **vars; // array of pointers to those variables: this must be allocated by user

    // Pointer to input-data stored in Chapter struct. For
    // documentation see the section.h file in src/etc directory. This
    // data-structure has all data read in from the input file and can
    // be used by applications to set application specific variables.
    Chapter *inpdata;

    // following data variables SHOULD NOT be set/modified by the
    // user: they are used internally

    double dx; // grid spacing
    int ncoeffs; // no of basis functions
    double nv[10]; // used to communicate some values

    int p_nipar;
    int *p_ipar; // algorithm specific internal integer variables
    
    int p_nrpar;
    double *p_rpar; // algorithm specific internal real variables

    FArray<double> kappa; // capacity function

    // pointer to workspace: this points to object of type
    // WAVE_Workspace or RKDG_Workspace
    void *work;

    // function pointer for core stepping function
    stepf_t step;
};

//
// function prototypes for standard functions
//
void solve(Run_Data& rd);
void advance(Run_Data& rd, FArray<double>& q, double tcurr, double tend);
void bc(const Run_Data& rd, FArray<double>& q);
void out(const Run_Data& rd, char const *fname, const FArray<double>& q);
int read_inp(Run_Data& rd, char *fname);
void driver(Run_Data& rd, char *fname);
void setprob(Run_Data& rd);
int init_output(Run_Data& ri);
void eval_fluctuations(Run_Data& rd, const FArray<double>& wave, const FArray<double>& s,
                       FArray<double>& amdq, FArray<double>& apdq);
void eval_fluctuations_qwave(Run_Data& rd, const FArray<double>& wave, const FArray<double>& s,
                             FArray<double>& amdq, FArray<double>& apdq);
void eval_fluctuations_fwave(Run_Data& rd, const FArray<double>& wave, const FArray<double>& s,
                             FArray<double>& amdq, FArray<double>& apdq);
void grid_transform(const Run_Data& rd, const FArray<double>& xc, FArray<double>& xp);
void write_grid(const Run_Data& rd, const FArray<double>& xp);
void before_step(Run_Data& rd, FArray<double>& q, double t);
void after_step(Run_Data& rd, FArray<double>& q, double t, double dt);

//
// prototypes for WAVE algorithm
//
void wave_setup(Run_Data& rd);
void wave_initialize(Run_Data& rd, FArray<double>& q);
void wave_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla);
void wave_limiter(Run_Data& rd, FArray<double>& wave, FArray<double>& s);

//
// prototypes for RKGD algorithms
//
void rkdg_setup(Run_Data& rd);
void rkdg_initialize(Run_Data& rd, FArray<double>& q);
void rkdg_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla);
void rkdg_limiter(Run_Data& rd, FArray<double>& q, double dt);

//
// prototypes for RKGD algorithms
//
void rkdg_setup(Run_Data& rd);
void rkdg_initialize(Run_Data& rd, FArray<double>& q);
void rkdg_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla);
void rkdg_limiter(Run_Data& rd, FArray<double>& q, double dt);

//
// prototypes for MACCOR2 algorithms
//
void maccor2_setup(Run_Data& rd);
void maccor2_initialize(Run_Data& rd, FArray<double>& q);
void maccor2_step(Run_Data& rd, FArray<double>& q, double tcurr, double dt, double& cfla);

#endif // __miniwarpx__2005__
