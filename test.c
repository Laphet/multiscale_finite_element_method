#include "math.h"
#include "sfem.h"

double a1(double x, double y) { return 1.0; }
double a2(double x, double y) { return 1.0; }
double a3(double x, double y) { return 0.0; }

double f(double x, double y) { return 2.0; }

double bdry(double x, double y) { return -1.0 * x * x; }

double u(double x, double y) { return -1.0 * x * x; }

int main(void)
{
    coefficient A = { .a1 = a1, .a2 = a2, .a3 = a3 };
    sfemInit(10);
    setCoefficient(A);
    int result = solvePDE(f, bdry);
    printf("error=%f\n", getError(u));
    sfemFinal();
    return result;
}
