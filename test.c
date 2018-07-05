#include "sfem.h"

double a1(double x, double y) { return 1.0; }
double a2(double x, double y) { return 2.0; }
double a3(double x, double y) { return 0.0; }

double f(double x, double y) { return 1.0; }

double test2(double x, double y) { return sin(x * y / 0.01); }

double test3(double x, double y) { return sqrt(x * y / 0.1); }

int main(void) 
{
    coefficient A = {.a1=a1, .a2=a2, .a3=a3};
    sfemInit(3);
    setCoefficient(A);
    int result = solveEPDE(f, zero);
    
    return 0; 
}
