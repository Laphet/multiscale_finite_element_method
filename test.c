#include "sfem.h"

double a1(double x, double y) { return 1.0; }
double a2(double x, double y) { return 2.0; }
double a3(double x, double y) { return 0.0; }

double test1(double x, double y)
{
    return cos((x+y+0.1)/0.1); 
}

double test2(double x, double y)
{
    return sin(x*y/0.1);
}

double test3(double x, double y)
{
    return sqrt(x*y/0.1);
}

int main(void)
{
    rectangle refRec = {.xmin = -1.0, .xmax = 1.0, .ymin = -1.0, .ymax = 1.0};
    printf("test1=%f\n", getNumericalIntegration(test1, refRec));
    return 0;
}
