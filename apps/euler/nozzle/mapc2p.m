clear all;

mx = 100;
x0 = 4.8;     
deg = 0;      

xl = 0; 
xu = 10;

a    = 10.0**-deg;
num1 = (xu-x0)/a;
num2 = (xl-x0)/a;
R    = log(num1+sqrt(num1**2+1.0))/ log(num2+sqrt(num2**2+1.0));
i0   = (mx-R)/(1.0-R);
b    = log(num2+sqrt(num2**2+1.0))/(1-i0);
c    = -b*i0;
d    = x0;

xp(1:100) = a*sinh(b*(1:100)+c)+d;
