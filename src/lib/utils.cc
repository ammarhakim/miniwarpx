#include <iostream>
#include <fstream>
#include <math.h>

#include "farray.h"
#include "utils.h"

double 
dmax(double a, double b)
{ return a>b?a:b; }

double 
dmax(double a, double b, double c)
{ 
    double m1= a>b?a:b;
    return m1>c?m1:c;
}

double 
dmin(double a, double b)
{ return a<b?a:b; }

double 
dmin(double a, double b, double c)
{ 
    double m1= a<b?a:b;
    return m1<c?m1:c;
}

double 
dsign(double a, double b)
{
    double sgn = (b>=0.0)?1.0:-1.0;
    return fabs(a)*sgn;
}

double 
mminmod(double a, double b, double c, double dx, double M)
{
    if (fabs(a) < M*dx*dx)
        return a;
    double sa = SGN(a);
    double sb = SGN(b);
    double sc = SGN(c);
    if( (sa==sb) && (sb==sc) )
    {
        if (sa<0)
            return dmax(a,b,c);
        else
            return dmin(a,b,c);
    }
    else
        return 0;
}

