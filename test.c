#include "math.h"
#include "sfem.h"

double a1(double x, double y) { return 1.0; }
double a2(double x, double y) { return 1.0; }
double a3(double x, double y) { return 0.5; }

double f(double x, double y) { return 2.0; }

double f1(double x, double y) { return -x; }

double f2(double x, double y) { return -y; }

double bdry(double x, double y) { return -2.0 * x * y; }

double u(double x, double y) { return -2.0 * x * y; }

int main(void)
{
    coefficient A = { .a1 = a1, .a2 = a2, .a3 = a3 };
    sfemInit(3);
    setCoefficient(A);
    int result = solvePDEwithDivF(f1, f2, bdry);
    sfemTest();
    printf("error=%f\n", getError(u));
    result = solvePDE(f, bdry);
    sfemTest();
    printf("error=%f\n", getError(u));
    sfemFinal();
    return result;
}
