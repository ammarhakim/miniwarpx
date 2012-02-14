#include <stdio.h>

double
stress(double bet, double modulus, double strain)
{
    double stres;

    stres = modulus*strain + bet*modulus*modulus*strain*strain;
    return stres;
}
