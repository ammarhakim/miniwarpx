#ifndef __rkdg_algo__2005__
#define __rkdg_algo__2005__

#include "miniwarpx.h"

//
// This file is part of WarpX hyperbolic conservation law library.
//
// Author: Ammar Hakim
// Date  : November 2005
//

/*
  Data and functions specific to the RKDG algorithm
 */

// the following transforms are coded as macros as they are called
// many many times

// ETA maps x -> eta
// ATE maps eta -> x

#define ETA(x,xcell,dx) (2.0*((x)-(xcell))/(dx))
#define ATE(eta,xcell,dx)((eta)/2.0*(dx)+(xcell))

/*
  Weights, w[0..n-1], and abscissa, x[0..n-1], for nth order Gaussian quadrature in
  interval x1..x2
*/
void gauleg(int n, double x1, double x2, FArray<double>& x, FArray<double>& w);

/**
   Legendre Polynomials for 'x' in [-1,1] and for given order 'n'.
 */
double legendre_poly(int n, double x);

/**
   Derivative of the Legendre Polynomials for 'x' in [-1,1] and for
   given order 'n'
 */
double legendre_poly_deriv(int n, double x);

void rkdg_eval_expansion(const Run_Data& rd, FArray<double>& qe, const FArray<double>& qi, int ec);
void rkdg_eval_expansion_left_edge(const Run_Data& rd, FArray<double>& qe, const FArray<double>& qi);
void rkdg_eval_expansion_right_edge(const Run_Data& rd, FArray<double>& qe, const FArray<double>& qi);

#endif // __rkdg_algo__2005__
