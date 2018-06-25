#include <math.h>
#include "gsl/gsl_integration.h"

typedef 
struct rectangle
{
    double xmin;
    double xmax;
    double ymin;
    double ymax;
};

typedef double (*funcP) (double, double);


const int LIMIT = 1000;
const double EPSABS = 0.0;
const double EPSREL = 1e-7;

funcP g_funcP = NULL;
gsl_integration_workspace *g_workspace = NULL;
rectangle g_rectangle = {.xmin = 0.0, .xmax = 1.0, .ymin = 0.0, .ymax = 1.0};


double getNumericalIntegration (funcP f, rectangle rectangleDom);



