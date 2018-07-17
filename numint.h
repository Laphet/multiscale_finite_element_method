#ifndef NUMINT_H
#define NUMINT_H

#include "gsl/gsl_integration.h"
#include "math.h"
#include "stdio.h"

typedef struct rectangle {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
} rectangle;
extern const int LIMIT;
extern const int NUM_POINT;
extern const double EPSABS;

extern const rectangle refRec;

typedef double (*func)(double, double);

double getNumericalIntegration(func f, rectangle rectangleDom);

#endif
