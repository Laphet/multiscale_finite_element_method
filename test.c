#include <stdio.h>
#include "numint.h"


double
f (double x, double y)
{
    return x * x * y;
}

int
main (void)
{
    rectangle rect = {.xmin = 2.0, .xmax = 3.0, .ymin = 2.0, .ymax = 3.0};
    printf ("the integration is %f.\n", getNumericalIntegration (f, rect));
    return 0;
}
