#include "math.h"
#include "sfem.h"

double a1(double x, double y) { return 1.0; }
double a2(double x, double y) { return 1.0; }
double a3(double x, double y) { return 0.5; }

double f(double x, double y) { return -1.0; }

double f1(double x, double y) { return -x; }

double f2(double x, double y) { return -y; }

double bdry(double x, double y) { return x*y; }

double u(double x, double y) { return x * y; }

int main(void)
{
    coefficient A = { .a1 = a1, .a2 = a2, .a3 = a3 };
    sfemInit(3);
    int result = solvePDE(A, f, bdry);
    printf("result is %d\n", result);
    sfemTest();
    printf("error is %.6f\n", getError(u));
    sfemFinal();
    return result;
}
