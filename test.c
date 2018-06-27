#include "sfem.h"

double a1(double x, double y) { return 1.0; }
double a2(double x, double y) { return 2.0; }
double a3(double x, double y) { return 0.0; }

int main(void)
{
    int N = 3;
    coefficient A = {.a1 = a1, .a2 = a2, .a3 = a3};
    sfemInit(N);
    setCoefficient(A);
    sfemFinal();
    return 0;
}
