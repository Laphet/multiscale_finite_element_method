#include "numint.h"

const int LIMIT = 10;
const int NUM_POINT = 10;
const double EPSABS = 1e-7;

func g_func = NULL;
gsl_integration_glfixed_table* g_glfix = NULL;
rectangle g_rectangle = { .xmin = 0.0, .xmax = 1.0, .ymin = 0.0, .ymax = 1.0 };

double gslInnerIntegrand(double y, void* p)
{
    double x = *(double*)p;
    return g_func(x, y);
}

double gslOuterIntergrand(double x, void* p)
{
    gsl_function F;
    F.function = &gslInnerIntegrand;
    F.params = &x;
    double resultPrev = 0.0, resultUpdated = 0.0;
    int iter = 0, i = 0, intervalCount = 1;
    double deltaY = 0.0, intervalStart = 0.0;
    while (1) {
        deltaY = (g_rectangle.ymax - g_rectangle.ymin) / intervalCount;
        intervalStart = g_rectangle.ymin;
        resultPrev = resultUpdated;
        resultUpdated = 0.0;
        for (i = 0; i < intervalCount; i++) {
            resultUpdated += gsl_integration_glfixed(
                &F, intervalStart, intervalStart + deltaY, g_glfix);
            intervalStart += deltaY;
        }
        if (fabs(resultPrev - resultUpdated) < EPSABS && iter < LIMIT)
            break;
        else if (iter >= LIMIT) {
            printf("reach iteration LIMIT(%d) \n at x=%f", LIMIT, x);
            break;
        } else {
            intervalCount = intervalCount * 2;
            iter++;
        }
    }
    return resultUpdated;
}

double getNumericalIntegration(func f, rectangle rectangleDom)
{
    g_func = f;
    g_rectangle = rectangleDom;
    g_glfix = gsl_integration_glfixed_table_alloc(NUM_POINT);
    gsl_function F;
    F.function = &gslOuterIntergrand;
    F.params = NULL;
    double resultPrev = 0.0, resultUpdated = 0.0;
    int iter = 0, i = 0, intervalCount = 1;
    double deltaX = 0.0, intervalStart = 0.0;
    while (1) {
        deltaX = (g_rectangle.xmax - g_rectangle.xmin) / intervalCount;
        intervalStart = g_rectangle.xmin;
        resultPrev = resultUpdated;
        resultUpdated = 0.0;
        for (i = 0; i < intervalCount; i++) {
            resultUpdated += gsl_integration_glfixed(
                &F, intervalStart, intervalStart + deltaX, g_glfix);
            intervalStart += deltaX;
        }
        if (fabs(resultPrev - resultUpdated) < EPSABS && iter < LIMIT)
            break;
        else if (iter >= LIMIT) {
            printf("reach iteration LIMIT(%d) \n", LIMIT);
            break;
        } else {
            intervalCount = intervalCount * 2;
            iter++;
        }
    }
    gsl_integration_glfixed_table_free(g_glfix);
    return resultUpdated;
}
