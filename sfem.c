#include "sfem.h"

int g_sliceNum = 2;
coefficient g_A = NULL;
double ***g_inteValueOnElement = NULL;
func g_f = NULL, g_f1 = NULL, g_f2 = NULL, g_bdry = NULL;
gsl_matrix *g_nodeValue = NULL;
gsl_matrix *g_stiffnessMatrix = NULL;
gsl_matrix *g_innerNodeValue = NULL;
gsl_vector *g_loadVector = NULL;
int g_m = 0, g_n = 0, g_flag = 0;
rectangle refRec = {.xmin = -1.0, .xmax = 1.0, .ymin = -1.0, .ymax = 1.0};
const int NUM_OF_INTE = 10;

int getIndexInVector(const int m, const int n)
{
    return (m - 1) * (g_sliceNum - 1) + n - 1;
}

double integrand(const double x, const double y)
{
    double h = 1.0 / g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag)
    {
    case 0:
        result = (y * y - 1.0) * g_A.a1(x_, y_) + (2 * x * y - 2.0) * g_A.a3(x_, y_) + (x * x - 1.0) * g_A.a2(x_, y_);
    case 1:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_) + 2 * (y - 1.0) * (x - 1.0) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
    case 2:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_) + 2 * (y - 1.0) * (x + 1.0) * g_A.a3(x_, y_) + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
    case 3:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_) + 2 * (y + 1.0) * (x + 1.0) * g_A.a3(x_, y_) + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
    case 4:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_) + 2 * (y + 1.0) * (x - 1.0) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
    case 5:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_) + ((x - 1.0) * (y - 1.0) + (y - 1.0) * (x + 1.0)) * g_A.a3(x_, y_) + (x - 1.0) * (x + 1.0) * g_A.a2(x_, y_);
    case 6:
        result = (y - 1.0) * (y + 1.0) * g_A.a1(x_, y_) + ((x + 1.0) * (y + 1.0) + (y - 1.0) * (x + 1.0)) * g_A.a3(x_, y_) + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
    case 7:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_) + ((x + 1.0) * (y + 1.0) + (y + 1.0) * (x - 1.0)) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
    case 8:
        result = (y + 1.0) * (y - 1.0) * g_A.a1(x_, y_) + ((x - 1.0) * (y - 1.0) + (y + 1.0) * (x - 1.0)) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
    case 9:
        result = (y * y - 1.0) * g_A.a1(x_, y_) + (2 * x * y + 2.0) * g_A.a3(x_, y_) + (x * x - 1.0) * g_A.a2(x_, y_);
    default:
        result = 0.0;
    }
    return h * h * result / 16.0;
}

func getIntegrand(int m, int n, int flag)
{
    g_m = m;
    g_n = n;
    g_flag = flag;
    return integrand;
}

void sfemInit(int sliceNum)
{
    g_sliceNum = sliceNum;
    ***g_inteValueOnElement = (double ***)malloc(sizeof(**double) * g_sliceNum);
    for (int i = 0; i < g_sliceNum; i++)
    {
        g_inteValueOnElement[i] = (double **)malloc(sizeof(*double) * g_sliceNum);
        for (int j = 0; j < g_sliceNum; j++)
            g_inteValueOnElement[i][j] = (double *)malloc(sizeof(double) * NUM_OF_INTE);
    }
    g_nodeValue = gsl_matrix_calloc(g_sliceNum + 1, g_sliceNum + 1);
    g_stiffnessMatrix = gsl_spmatrix_alloc ((g_sliceNum-1) * (g_sliceNum-1), (g_sliceNum-1) * (g_sliceNum-1);
    g_innerNodeValue = gsl_vector_calloc ((g_sliceNum-1) * (g_sliceNum-1));
    g_loadVector = gsl_vector_calloc ((g_sliceNum-1) * (g_sliceNum-1));
}

void setCoefficient(coefficient A)
{
    g_A = A;
    for (int m = 0; m < g_sliceNum; m++)
        for (int n = 0; n < g_sliceNum; n++)
            for (int flag = 0; flag < NUM_OF_INTE; flag++)
                g_inteValueOnElement[m][n][flag] = getNumericalIntegration(getIntegrand(m, n, flag), refRec);

    /* test */
    for (int m = 0; m < g_sliceNum; m++)
        for (int n = 0; n < g_sliceNum; n++)
        {
            printf('\n m=%d n=%d :', m, n);
            for (int flag = 0; flag < NUM_OF_INTE; flag++)
                printf('%f\t', g_inteValueOnElement[m][n][flag]);
        }
}

void solveEPDE(func f, func bdry)
{
    /*
    g_A = A;
    g_f = f;
    g_bdry = bdry;
    for (int m = 0; m < g_sliceNum + 1; m++)
    {
        for (int n = m; n < g_sliceNum + 1; n++)
        {
            double x = (double) m / g_sliceNum;
            double y = (double) n / g_sliceNum;
            if (m == 0 || n == 0 || m == g_sliceNum || n == g_sliceNum)
                gsl_matrix_set (g_nodeValue, m, n, g_bdry (x, y));
            elseif (m == 1 || n == 1 || m == g_sliceNum - 1 || n == g_sliceNum - 1)
            {


            }
            else
            {

            }




        }
    }
    */
}

void solveEPDEwithDivF(func f1, func f2, func bdry)
{
}

void sfemFinal(void)
{
    gsl_block_free(g_nodeValue);
    gsl_spmatrix_free(g_stiffnessMatrix);
    gsl_vector_free(g_innerNodeValue);
    gsl_vector_free(g_loadVector);
    for (int i = 0; i < g_sliceNum; i++)
    {
        for (int j = 0; j < g_sliceNum; j++)
            free(g_inteValueOnElement[i][j]);
        free(g_inteValueOnElement[i]);
    }
    free(g_inteValueOnElement);
}
