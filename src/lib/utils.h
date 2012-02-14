#ifndef __utils_2005__
#define __utils_2005__

#ifndef PI
#define PI 3.14159265358979323846 
#endif

#include <math.h>

#define SQR(x) pow((x),2)
#define SGN(b) ((b)>=0.)?1.:-1.

/*
 * Utility functions
 */
// max of two doubles
double dmax(double a, double b);
// max of three doubles
double dmax(double a, double b, double c);

// min of two doubles
double dmin(double a, double b);
// min of three doubles
double dmin(double a, double b, double c);

// |a|*sign(b)
double dsign(double a, double b);

// Modified min-mod function
double mminmod(double a, double b, double c, double dx, double M);

#endif // __utils_2005__
