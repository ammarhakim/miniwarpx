#include "rkdg_algo.h"
#include "utils.h"
#include <math.h>

/**
  Weights, w, and abscissa, x, for nth order Gaussian quadrature in
  interval x1..x2.
 */
void 
gauleg(int n, double x1, double x2, FArray<double>& x, FArray<double>& w)
{
    // call routine adapted from Numerical Recipies in C book
    extern void __gauleg( double x1, double x2,  double x[], double w[], int n);
    // note that as __gauleg expects ordinary C arrays indexed 1..n we
    // need to subtract 1 from data()
    __gauleg(x1, x2, x.data()-1, w.data()-1, n);
}


/**
   Legendre Polynomials for 'x' in [-1,1] and for given order 'n'.
 */
double
legendre_poly(int n, double x)
{
    double p0 = 1.0;
    double p1 = x;
    double pn,pn1,pn2;

    if (n==0) return p0;
    if (n==1) return p1;

    pn = 0.0; pn1 = p1; pn2 = p0; // initialize recurrence
    for(int i=2; i<=n; i++)
    {
        // use recurrence relation to compute P_n
        pn = (x*(2.*i-1.)*pn1 - (i-1.)*pn2)/(1.*i);
        pn2 = pn1;
        pn1 = pn;
    }
    return pn;
}

/**
   Derivative of the Legendre Polynomials for 'x' in [-1,1] and for
   given order 'n'
 */
double 
legendre_poly_deriv(int n, double x)
{
    double dp0 = 0;
    double dp1 = 1;
    double dpn;

    if (n==0) return dp0;
    if (n==1) return dp1;

    dpn = ((n+1.)*x*legendre_poly(n,x) - (n+1.)*legendre_poly(n+1,x))/(1.-x*x);

    return dpn;
}
