#include <stdio.h>
#include <math.h>
#include "gsl/gsl_integration.h"


const int LIMIT = 1000;
double g_epsabs = 0.0;
double g_epsrel = 1e-7;
double (*g_integrandFuncP) (double, double) = NULL;
gsl_integration_workspace * g_workspace = NULL;
double g_xmin = 0.0, g_xmax = 0.0, g_ymin = 0.0, g_ymax = 0.0;


double
gslInnerIntegrand (double y, void * p)
{
    double x = *(double *) p;
    return g_integrandFuncP (x, y);
}


double
gslOuterIntergrand (double x, void * p)
{
    gsl_function F;
    F.function = &gslInnerIntegrand;
    F.params = & x;
    double result = 0.0, error = 0.0;
    gsl_integration_qag (&F, g_ymin, g_ymax, g_epsabs, g_epsrel, LIMIT, 1, g_workspace, &result, &error);
    return result;
}


double
getNumericalIntegration (const double (*integrandFuncP) (double, double), const double xmin, const double xmax,
                            const double ymin, const double ymax)
{
    g_integrandFuncP = integrandFuncP;
    g_xmin = xmin;
    g_xmax = xmax;
    g_ymin = ymin;
    g_ymax = ymax;
    g_workspace = gsl_integration_workspace_alloc (LIMIT);
    gsl_function F;
    F.function = &gslOuterIntergrand;
    F.params = NULL;
    double result = 0.0, error = 0.0;
    gsl_integration_qag (&F, g_xmin, g_xmax, g_epsabs, g_epsrel, LIMIT, 1, g_workspace, &result, &error);
    gsl_integration_workspace_free (g_workspace);
    return result;
}


double
integrand (double x, double y)
{
    return 1.0
}


int 
main(void)
{
    double xmin = 0.0, xmax = 1.0, ymin  = 0.0, ymax = 1.0;
    double result = getNumericalIntegration (&integrand, xmin, xmax, ymin, ymax);
    printf ("the Integral Value = %f.\n", result);
    return 0;

}
