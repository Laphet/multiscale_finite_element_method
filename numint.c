#include "numint.h"


double
gslInnerIntegrand (double y, void * p)
{
    double x = *(double *) p;
    return g_funcP (x, y);
}


double
gslOuterIntergrand (double x, void * p)
{
    gsl_function F;
    F.function = &gslInnerIntegrand;
    F.params = &x;
    double result = 0.0, error = 0.0;
    gsl_integration_qag (&F, g_rectangle.ymin, grectangle.ymax, EPSABS, EPSREL, LIMIT, 1, g_workspace, &result, &error);
    return result;
}


double
getNumericalIntegration (funcP f, rectangle rectangleDom)
{
    g_funcP = f;
    g_rectangle = rectangleDom;
    g_workspace = gsl_integration_workspace_alloc (LIMIT);
    gsl_function F;
    F.function = &gslOuterIntergrand;
    F.params = NULL;
    double result = 0.0, error = 0.0;
    gsl_integration_qag (&F, g_rectangle.xmin, g_rectangle.xmax, EPSABS, EPSREL, LIMIT, 1, g_workspace, &result, &error);
    gsl_integration_workspace_free (g_workspace);
    return result;
}