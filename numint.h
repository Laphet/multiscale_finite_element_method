#ifndef NUMINT_H
#define NUMINT_H

#include "gsl/gsl_integration.h"

typedef 
struct rectangle
{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
}rectangle;

typedef double (*funcP) (double, double);

double getNumericalIntegration (funcP f, rectangle rectangleDom);

#endif

