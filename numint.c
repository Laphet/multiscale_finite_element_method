#include "numint.h"

const int LIMIT = 1000;
const double EPSABS = 0.0;
const double EPSREL = 1e-7;

func g_func = NULL;
gsl_integration_workspace *g_workspace = NULL;
rectangle g_rectangle = {.xmin = 0.0, .xmax = 1.0, .ymin = 0.0, .ymax = 1.0};

double
gslInnerIntegrand(double y, void *p)
{
    double x = *(double *)p;
    return g_func(x, y);
}

double
gslOuterIntergrand(double x, void *p)
{
    gsl_function F;
    F.function = &gslInnerIntegrand;
    F.params = &x;
    double result = 0.0, error = 0.0;
    gsl_integration_qag(&F, g_rectangle.ymin, g_rectangle.ymax, EPSABS, EPSREL, LIMIT, 1, g_workspace, &result, &error);
    return result;
}

double
getNumericalIntegration(func f, rectangle rectangleDom)
{
    g_func = f;
    g_rectangle = rectangleDom;
    g_workspace = gsl_integration_workspace_alloc(LIMIT);
    gsl_function F;
    F.function = &gslOuterIntergrand;
    F.params = NULL;
    double result = 0.0, error = 0.0;
    gsl_integration_qag(&F, g_rectangle.xmin, g_rectangle.xmax, EPSABS, EPSREL, LIMIT, 1, g_workspace, &result, &error);
    gsl_integration_workspace_free(g_workspace);
    return result;
}
