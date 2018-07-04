#include "sfem.h"

int g_sliceNum = 2;
coefficient g_A;
double ***g_inteValueOnElement = NULL, ***g_loadValueOnElement = NULL;
func g_f = NULL, g_f1 = NULL, g_f2 = NULL, g_bdry = NULL;
gsl_matrix *g_nodeValue = NULL;
gsl_spmatrix *g_stiffnessMatrix = NULL;
gsl_spmatrix *g_stiffnessMatrixExt = NULL;
gsl_vector *g_innerNodeValue = NULL;
gsl_vector *g_loadVector = NULL;
int g_m = 0, g_n = 0, g_flag = 0;
rectangle refRec = {.xmin = -1.0, .xmax = 1.0, .ymin = -1.0, .ymax = 1.0};
const int NUM_OF_INTE = 10, NUM_OF_LOAD = 4;

double zero(double x, double y) { return 0.0; }

int getIndexInVector(const int m, const int n)
{
    return n * (g_sliceNum + 1) + m;
}

int *getIndexInMatrix(const int index)
{
    int *entry = (int *)malloc(sizeof(int) * 2);
    entry[0] = index % (g_sliceNum + 1);
    entry[1] = index / (g_sliceNum + 1);
    return entry;
}

double integrand(double x, double y)
{
    double h = 1.0 / (double)g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag)
    {
    case 0:
        result = (y * y - 1.0) * g_A.a1(x_, y_) + (2 * x * y - 2.0) * g_A.a3(x_, y_) + (x * x - 1.0) * g_A.a2(x_, y_);
        break;
    case 1:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_) + 2 * (y - 1.0) * (x - 1.0) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
        break;
    case 2:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_) + 2 * (y - 1.0) * (x + 1.0) * g_A.a3(x_, y_) + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
        break;
    case 3:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_) + 2 * (y + 1.0) * (x + 1.0) * g_A.a3(x_, y_) + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
        break;
    case 4:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_) + 2 * (y + 1.0) * (x - 1.0) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
        break;
    case 5:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_) + ((x - 1.0) * (y - 1.0) + (y - 1.0) * (x + 1.0)) * g_A.a3(x_, y_) + (x - 1.0) * (x + 1.0) * g_A.a2(x_, y_);
        break;
    case 6:
        result = (y - 1.0) * (y + 1.0) * g_A.a1(x_, y_) + ((x + 1.0) * (y + 1.0) + (y - 1.0) * (x + 1.0)) * g_A.a3(x_, y_) + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
        break;
    case 7:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_) + ((x + 1.0) * (y + 1.0) + (y + 1.0) * (x - 1.0)) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
        break;
    case 8:
        result = (y + 1.0) * (y - 1.0) * g_A.a1(x_, y_) + ((x - 1.0) * (y - 1.0) + (y + 1.0) * (x - 1.0)) * g_A.a3(x_, y_) + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
        break;
    case 9:
        result = (y * y - 1.0) * g_A.a1(x_, y_) + (2 * x * y + 2.0) * g_A.a3(x_, y_) + (x * x - 1.0) * g_A.a2(x_, y_);
        break;
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

double integrandOfLoad(double x, double y)
{
    double h = 1.0 / (double)g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag)
    {
    case 0:
        result = g_f(x_, y_) * (x - 1.0) * (y - 1.0);
        break;
    case 1:
        result = g_f(x_, y_) * (x + 1.0) * (y - 1.0);
        break;
    case 2:
        result = g_f(x_, y_) * (x + 1.0) * (y + 1.0);
        break;
    case 3:
        result = g_f(x_, y_) * (x - 1.0) * (y + 1.0);
    default:
        result = 0.0;
    }
    return h * h * result / 4.0;
}

func getIntegrandOfLoad(int m, int n, int flag)
{
    g_m = m;
    g_n = n;
    g_flag = flag;
    return integrandOfLoad;
}

double integrandOfLoadWithDivF(double x, double y)
{
    double h = 1.0 / (double)g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag)
    {
    case 0:
        result = g_f(x_, y_) * (x - 1.0) * (y - 1.0);
        break;
    case 1:
        result = g_f(x_, y_) * (x + 1.0) * (y - 1.0);
        break;
    case 2:
        result = g_f(x_, y_) * (x + 1.0) * (y + 1.0);
        break;
    case 3:
        result = g_f(x_, y_) * (x - 1.0) * (y + 1.0);
    default:
        result = 0.0;
    }
    return h * h * result / 4.0;
}

func getIntegrandOfLoadWithDivF(int m, int n, int flag)
{
    g_m = m;
    g_n = n;
    g_flag = flag;
    return integrandOfLoadWithDivF;
}

void setLoadVector(int isWithDivF)
{
    int m = 0, n = 0, flag = 0, i = 0;
    int *entry = NULL;
    double temp = 0.0;
    for (m = 0; m < g_sliceNum; m++)
        for (n = 0; n < g_sliceNum; n++)
            for (flag = 0; flag < NUM_OF_LOAD, flag++)
                if (isWithDivF)
                    g_loadValueOnElement[m][n][flag] = getNumericalIntegration(getIntegrandOfLoadWithDivF(m, n, flag), refRec);
                else
                    g_loadValueOnElement[m][n][flag] = getNumericalIntegration(getIntegrandOfLoad(m, n, flag), refRec);

    for (i = 0; i < (g_sliceNum - 1) * (g_sliceNum - 1); i++)
    {
        m = i % (g_sliceNum - 1) + 1;
        n = i / (g_sliceNum - 1) + 1;
        temp = g_loadValueOnElement[m][n][0] + g_loadValueOnElement[m - 1][n][1] +
               g_loadValueOnElement[m - 1][n - 1][2] + g_loadValueOnElement[m][n - 1][3];
        gsl_vector_set(g_loadVector, i, temp);
    }
}

void containerInit(int num, double ***container)
{
    int m = 0, n = 0;
    container = (double ***)malloc(sizeof(double **) * g_sliceNum);
    for (m = 0; m < g_sliceNum; m++)
    {
        container[m] = (double **)malloc(sizeof(double *) * g_sliceNum);
        for (n = 0; n < g_sliceNum; n++)
            container[m][n] = (double *)malloc(sizeof(double) * num);
    }
}

void containerFree(double ***container)
{
    int m = 0, n = 0;
    for (m = 0; m < g_sliceNum; m++)
    {
        for (n = 0; n < g_sliceNum; n++)
            free(container[m][n]);
        free(container[m]);
    }
    free(container);
}

void sfemInit(int sliceNum)
{
    g_sliceNum = sliceNum;
    containerInit(NUM_OF_INTE, g_inteValueOnElement);
    containerInit(NUM_OF_LOAD, g_loadValueOnElement);
    g_nodeValue = gsl_matrix_calloc(g_sliceNum + 1, g_sliceNum + 1);
    g_stiffnessMatrix = gsl_spmatrix_alloc((g_sliceNum - 1) * (g_sliceNum - 1), (g_sliceNum - 1) * (g_sliceNum - 1));
    g_stiffnessMatrixExt = gsl_spmatrix_alloc((g_sliceNum + 1) * (g_sliceNum + 1), (g_sliceNum + 1) * (g_sliceNum + 1));
    g_innerNodeValue = gsl_vector_calloc((g_sliceNum - 1) * (g_sliceNum - 1));
    g_loadVector = gsl_vector_calloc((g_sliceNum - 1) * (g_sliceNum - 1));
}

void setCoefficient(coefficient A)
{
    g_A = A;
    int m = 0, n = 0, flag = 0;
    for (m = 0; m < g_sliceNum; m++)
        for (n = 0; n < g_sliceNum; n++)
            for (flag = 0; flag < NUM_OF_INTE; flag++)
                g_inteValueOnElement[m][n][flag] = getNumericalIntegration(getIntegrand(m, n, flag), refRec);
}

void setStiffnessMatrixExt(void)
{
    int i = 0, m = 0, n = 0, isInteriorPoint = 0;
    int *entry = NULL;
    double temp = 0.0;
    for (i = 0; i < (g_sliceNum + 1) * (g_sliceNum + 1); i++)
    {
        entry = getIndexInMatrix(i);
        m = entry[0];
        n = entry[1];
        isInteriorPoint = (m > 0) && (m < g_sliceNum) && (n > 0) && (n < g_sliceNum);
        if (isInteriorPoint)
        {
            temp = g_inteValueOnElement[m][n][1] + g_inteValueOnElement[m - 1][n][2] +
                   g_inteValueOnElement[m - 1][n - 1][3] + g_inteValueOnElement[m][n - 1][4];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m, n), temp);
            temp = g_inteValueOnElement[m][n - 1][7] + g_inteValueOnElement[m][n][5];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n), temp);
            temp = g_inteValueOnElement[m][n][8] + g_inteValueOnElement[m - 1][n][6];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m, n + 1), temp);
            temp = g_inteValueOnElement[m - 1][n][5] + g_inteValueOnElement[m - 1][n - 1][7];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m - 1, n), temp);
            temp = g_inteValueOnElement[m - 1][n - 1][6] + g_inteValueOnElement[m][n - 1][8];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m, n - 1), temp);
            temp = g_inteValueOnElement[m][n][0];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n + 1), temp);
            temp = g_inteValueOnElement[m - 1][n][9];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m - 1, n + 1), temp);
            temp = g_inteValueOnElement[m - 1][n - 1][0];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m - 1, n - 1), temp);
            temp = g_inteValueOnElement[m][n - 1][9];
            gsl_spmatrix_set(g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n - 1), temp);
        }
    }
}

void setLinearEquations2Solve()
{
    int i = 0, j = 0, m = 0, n = 0, isContactWithBdry = 0, b = 0;
    double temp = 0.0, h = 1.0 / (double)g_sliceNum;
    for (i = 0; i < (g_sliceNum - 1) * (g_sliceNum - 1); i++)
    {
        m = i % (g_sliceNum - 1) + 1;
        n = i / (g_sliceNum - 1) + 1;
        isContactWithBdry = (m == 1) || (m == g_sliceNum - 1) || (n == 1) || (n == g_sliceNum - 1);
        temp = 0.0;
        if (isContactWithBdry)
        {
            for (b = 1; b < g_sliceNum; b++)
            {
                temp += (g_bdry(0.0, (double)b * h) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(0, b)) +
                         g_bdry(1.0, (double)b * h) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(g_sliceNum, b)) +
                         g_bdry((double)b * h, 0.0) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(b, 0)) +
                         g_bdry((double)b * h, 1.0) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(b, g_sliceNum)));
            }
            temp += (g_bdry(0.0, 0.0) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(0, 0)) +
                     g_bdry(1.0, 0.0) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(g_sliceNum, 0)) +
                     g_bdry(1.0, 1.0) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(g_sliceNum, g_sliceNum)) +
                     g_bdry(0.0, 1.0) * gsl_spmatrix_get(g_stiffnessMatrixExt, i, getIndexInVector(0, g_sliceNum)));
            gsl_vector_set(g_loadVector, i, gsl_vector_get(g_loadVector, i) - temp);
        }
        for (j = 0; j < i + 1; j++)
        {
            gsl_spmatrix_set(g_stiffnessMatrix, i, j, gsl_spmatrix_get(g_stiffnessMatrixExt, i + 1, j + 1));
            gsl_spmatrix_set(g_stiffnessMatrix, j, i, gsl_spmatrix_get(g_stiffnessMatrixExt, j + 1, i + 1));
        }
    }
}

int solveEPDE(func f, func bdry)
{
    g_f = f;
    g_bdry = bdry;
    setStiffnessMatrixExt();
    setLoadVector(0);
    setLinearEquations2Solve();
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, (g_sliceNum-1)*(g_sliceNum-1), 5);
    return gsl_splinalg_itersolve_iterate(g_stiffnessMatrix, g_loadVector, ERR_TOL, g_innerNodeValue, work);
}

void solveEPDEwithDivF(func f1, func f2, func bdry)
{
}

void sfemFinal(void)
{
    gsl_matrix_free(g_nodeValue);
    gsl_spmatrix_free(g_stiffnessMatrix);
    gsl_spmatrix_free(g_stiffnessMatrixExt);
    gsl_vector_free(g_innerNodeValue);
    gsl_vector_free(g_loadVector);
    containerFree(g_inteValueOnElement);
    containerFree(g_loadValueOnElement);
}
