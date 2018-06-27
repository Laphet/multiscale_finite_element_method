#ifndef NUMINT_H
#define NUMINT_H

#include "gsl/gsl_integration.h"

typedef struct rectangle
{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
} rectangle;

typedef double (*func)(double, double);

double getNumericalIntegration(func f, rectangle rectangleDom);

#endif
