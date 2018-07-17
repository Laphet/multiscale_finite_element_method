#include "sfem.h"

static int g_sliceNum;
static coefficient g_A;
static func g_f, g_f1, g_f2, g_bdry;
static gsl_spmatrix* g_stiffnessMatrix;
static gsl_spmatrix* g_stiffnessMatrixExt;
static gsl_vector* g_innerNodeValue;
static gsl_vector* g_loadVector;
static int g_m, g_n, g_flag;

static const int NUM_OF_INTE = 10, NUM_OF_LOAD = 4;
static const double ERR_TOL = 1.0e-7;
static const int MAX_ITER = 1024;

void sfemInit(int sliceNum)
{
    g_sliceNum = sliceNum;
    g_stiffnessMatrix = gsl_spmatrix_alloc((g_sliceNum - 1) * (g_sliceNum - 1),
        (g_sliceNum - 1) * (g_sliceNum - 1));
    g_stiffnessMatrixExt
        = gsl_spmatrix_alloc((g_sliceNum - 1) * (g_sliceNum - 1),
            (g_sliceNum + 1) * (g_sliceNum + 1));
    g_innerNodeValue = gsl_vector_calloc((g_sliceNum - 1) * (g_sliceNum - 1));
    g_loadVector = gsl_vector_calloc((g_sliceNum - 1) * (g_sliceNum - 1));
}

static int getIndexInVector(const int m, const int n)
{
    return n * (g_sliceNum + 1) + m;
}

static double integrand(double x, double y)
{
    double h = 1.0 / (double)g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag) {
    case 0:
        result = (y * y - 1.0) * g_A.a1(x_, y_)
            + (2 * x * y - 2.0) * g_A.a3(x_, y_)
            + (x * x - 1.0) * g_A.a2(x_, y_);
        break;
    case 1:
        result = (y - 1.0) * (y - 1.0) * g_A.a1(x_, y_)
            + 2 * (y - 1.0) * (x - 1.0) * g_A.a3(x_, y_)
            + (x - 1.0) * (x - 1.0) * g_A.a2(x_, y_);
        break;
    case 2:
        result = (1.0 - y) * (1.0 - y) * g_A.a1(x_, y_)
            - 2 * (1.0 - y) * (1.0 + x) * g_A.a3(x_, y_)
            + (1.0 + x) * (1.0 + x) * g_A.a2(x_, y_);
        break;
    case 3:
        result = (y + 1.0) * (y + 1.0) * g_A.a1(x_, y_)
            + 2 * (y + 1.0) * (x + 1.0) * g_A.a3(x_, y_)
            + (x + 1.0) * (x + 1.0) * g_A.a2(x_, y_);
        break;
    case 4:
        result = (1.0 + y) * (1.0 + y) * g_A.a1(x_, y_)
            - 2 * (1.0 + y) * (1.0 - x) * g_A.a3(x_, y_)
            + (1.0 - x) * (1.0 - x) * g_A.a2(x_, y_);
        break;
    case 5:
        result = (y - 1.0) * (1.0 - y) * g_A.a1(x_, y_)
            + ((x - 1.0) * (1.0 - y) - (y - 1.0) * (1.0 + x)) * g_A.a3(x_, y_)
            - (x - 1.0) * (1.0 + x) * g_A.a2(x_, y_);
        break;
    case 6:
        result = (1.0 - y) * (y + 1.0) * g_A.a1(x_, y_)
            + (-(1.0 + x) * (y + 1.0) + (1.0 - y) * (x + 1.0)) * g_A.a3(x_, y_)
            - (1.0 + x) * (x + 1.0) * g_A.a2(x_, y_);
        break;
    case 7:
        result = -(y + 1.0) * (1.0 + y) * g_A.a1(x_, y_)
            + (-(x + 1.0) * (1.0 + y) + (y + 1.0) * (1.0 - x)) * g_A.a3(x_, y_)
            + (x + 1.0) * (1.0 - x) * g_A.a2(x_, y_);
        break;
    case 8:
        result = -(1.0 + y) * (y - 1.0) * g_A.a1(x_, y_)
            + ((1.0 - x) * (y - 1.0) - (1.0 + y) * (x - 1.0)) * g_A.a3(x_, y_)
            + (1.0 - x) * (x - 1.0) * g_A.a2(x_, y_);
        break;
    case 9:
        result = (y * y - 1.0) * g_A.a1(x_, y_)
            + (2 * x * y + 2.0) * g_A.a3(x_, y_)
            + (x * x - 1.0) * g_A.a2(x_, y_);
        break;
    default:
        result = 0.0;
    }
    return result / 16.0;
}

static double integrandOfLoadWithDivF(double x, double y)
{
    double h = 1.0 / (double)g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag) {
    case 0:
        result = g_f1(x_, y_) * (y - 1.0) + g_f2(x_, y_) * (x - 1.0);
        break;
    case 1:
        result = g_f1(x_, y_) * (1.0 - y) - g_f2(x_, y_) * (1.0 + x);
        break;
    case 2:
        result = g_f1(x_, y_) * (y + 1.0) + g_f2(x_, y_) * (x + 1.0);
        break;
    case 3:
        result = -g_f1(x_, y_) * (1.0 + y) + g_f2(x_, y_) * (1.0 - x);
        break;
    default:
        result = 0.0;
    }
    return result / 8.0;
}

static double integrandOfLoad(double x, double y)
{
    double h = 1.0 / (double)g_sliceNum;
    double x_ = ((x + 1.0) / 2.0 + (double)g_m) * h;
    double y_ = ((y + 1.0) / 2.0 + (double)g_n) * h;
    double result = 0.0;
    switch (g_flag) {
    case 0:
        result = g_f(x_, y_) * (1.0 - x) * (1.0 - y);
        break;
    case 1:
        result = g_f(x_, y_) * (1.0 + x) * (1.0 - y);
        break;
    case 2:
        result = g_f(x_, y_) * (1.0 + x) * (1.0 + y);
        break;
    case 3:
        result = g_f(x_, y_) * (1.0 - x) * (1.0 + y);
        break;
    default:
        result = 0.0;
    }
    return result / 16.0;
}

static func getIntegrand(int m, int n, int flag)
{
    g_m = m;
    g_n = n;
    g_flag = flag;
    return integrand;
}

static func getIntegrandOfLoad(int m, int n, int flag, int isWithDivF)
{
    g_m = m;
    g_n = n;
    g_flag = flag;
    if (isWithDivF)
        return integrandOfLoadWithDivF;
    else
        return integrandOfLoad;
}

static double getLocalStiffness(int m, int n, int flag)
{
    return getNumericalIntegration(getIntegrand(m, n, flag), refRec);
}

static double getLocalLoad(int m, int n, int flag, int isWithDivF)
{
    double h = 1.0 / (double)g_sliceNum;
    if (isWithDivF)
        return h
            * getNumericalIntegration(
                  getIntegrandOfLoad(m, n, flag, isWithDivF), refRec);
    else
        return h * h
            * getNumericalIntegration(
                  getIntegrandOfLoad(m, n, flag, isWithDivF), refRec);
}

static void setStiffnessMatrix(void)
{
    int i = 0, j = 0, m = 0, n = 0, m_ = 0, n_ = 0;
    double temp = 0.0;
    for (i = 0; i < (g_sliceNum - 1) * (g_sliceNum - 1); ++i) {
        m = i % (g_sliceNum - 1) + 1;
        n = i / (g_sliceNum - 1) + 1;
        temp = getLocalStiffness(m, n, 1) + getLocalStiffness(m - 1, n, 2)
            + getLocalStiffness(m - 1, n - 1, 3)
            + getLocalStiffness(m, n - 1, 4);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m, n), temp);
        temp = getLocalStiffness(m, n - 1, 7) + getLocalStiffness(m, n, 5);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n), temp);
        temp = getLocalStiffness(m, n, 8) + getLocalStiffness(m - 1, n, 6);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m, n + 1), temp);
        temp = getLocalStiffness(m - 1, n, 5)
            + getLocalStiffness(m - 1, n - 1, 7);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m - 1, n), temp);
        temp = getLocalStiffness(m - 1, n - 1, 6)
            + getLocalStiffness(m, n - 1, 8);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m, n - 1), temp);
        temp = getLocalStiffness(m, n, 0);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n + 1), temp);
        temp = getLocalStiffness(m - 1, n, 9);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m - 1, n + 1), temp);
        temp = getLocalStiffness(m - 1, n - 1, 0);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m - 1, n - 1), temp);
        temp = getLocalStiffness(m, n - 1, 9);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n - 1), temp);
        for (j = 0; j < i + 1; ++j) {
            m_ = j % (g_sliceNum - 1) + 1;
            n_ = j / (g_sliceNum - 1) + 1;
            temp = gsl_spmatrix_get(
                g_stiffnessMatrixExt, i, getIndexInVector(m_, n_));
            gsl_spmatrix_set(g_stiffnessMatrix, i, j, temp);
            gsl_spmatrix_set(g_stiffnessMatrix, j, i, temp);
        }
    }
}

static void setLoadVector(int isWithDivF)
{
    int m = 0, n = 0, i = 0, isContactWithBdry = 0, b = 0;
    double temp = 0.0, h = 1.0 / (double)g_sliceNum;
    for (i = 0; i < (g_sliceNum - 1) * (g_sliceNum - 1); ++i) {
        m = i % (g_sliceNum - 1) + 1;
        n = i / (g_sliceNum - 1) + 1;
        temp = getLocalLoad(m, n, 0, isWithDivF)
            + getLocalLoad(m - 1, n, 1, isWithDivF)
            + getLocalLoad(m - 1, n - 1, 2, isWithDivF)
            + getLocalLoad(m, n - 1, 3, isWithDivF);
        isContactWithBdry = (m == 1) || (m == g_sliceNum - 1) || (n == 1)
            || (n == g_sliceNum - 1);
        if (isContactWithBdry) {
            for (b = 1; b < g_sliceNum; ++b) {
                temp -= (g_bdry(0.0, (double)b * h)
                        * gsl_spmatrix_get(
                              g_stiffnessMatrixExt, i, getIndexInVector(0, b))
                    + g_bdry(1.0, (double)b * h)
                        * gsl_spmatrix_get(g_stiffnessMatrixExt, i,
                              getIndexInVector(g_sliceNum, b))
                    + g_bdry((double)b * h, 0.0)
                        * gsl_spmatrix_get(
                              g_stiffnessMatrixExt, i, getIndexInVector(b, 0))
                    + g_bdry((double)b * h, 1.0)
                        * gsl_spmatrix_get(g_stiffnessMatrixExt, i,
                              getIndexInVector(b, g_sliceNum)));
            }
            temp -= (g_bdry(0.0, 0.0)
                    * gsl_spmatrix_get(
                          g_stiffnessMatrixExt, i, getIndexInVector(0, 0))
                + g_bdry(1.0, 0.0)
                    * gsl_spmatrix_get(g_stiffnessMatrixExt, i,
                          getIndexInVector(g_sliceNum, 0))
                + g_bdry(1.0, 1.0)
                    * gsl_spmatrix_get(g_stiffnessMatrixExt, i,
                          getIndexInVector(g_sliceNum, g_sliceNum))
                + g_bdry(0.0, 1.0)
                    * gsl_spmatrix_get(g_stiffnessMatrixExt, i,
                          getIndexInVector(0, g_sliceNum)));
        }

        gsl_vector_set(g_loadVector, i, temp);
    }
}

int solvePDE(coefficient A, func f, func bdry)
{
    g_A = A;
    g_f = f;
    g_bdry = bdry;
    setStiffnessMatrix();
    setLoadVector(0);
    const gsl_splinalg_itersolve_type* T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve* work = gsl_splinalg_itersolve_alloc(
        T, (g_sliceNum - 1) * (g_sliceNum - 1), 0);
    int status = 0, iter = 0;
    do {
        status = gsl_splinalg_itersolve_iterate(
            g_stiffnessMatrix, g_loadVector, ERR_TOL, g_innerNodeValue, work);
        ++iter;
    } while (status == GSL_CONTINUE && iter < MAX_ITER);
    gsl_splinalg_itersolve_free(work);
    return 1;
}

int solvePDEwithDivF(coefficient A, func f1, func f2, func bdry)
{
    g_A = A;
    g_f1 = f1;
    g_f2 = f2;
    g_bdry = bdry;
    setStiffnessMatrix();
    setLoadVector(1);
    const gsl_splinalg_itersolve_type* T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve* work = gsl_splinalg_itersolve_alloc(
        T, (g_sliceNum - 1) * (g_sliceNum - 1), 0);
    int status = 0, iter = 0;
    do {
        status = gsl_splinalg_itersolve_iterate(
            g_stiffnessMatrix, g_loadVector, ERR_TOL, g_innerNodeValue, work);
        ++iter;
    } while (status == GSL_CONTINUE && iter < MAX_ITER);
    gsl_splinalg_itersolve_free(work);
    return 1;
}

gsl_vector* getInnerNodeValue()
{
    gsl_vector* innerNodeValue
        = gsl_vector_calloc((g_sliceNum - 1) * (g_sliceNum - 1));
    gsl_vector_memcpy(innerNodeValue, g_innerNodeValue);
    return innerNodeValue;
}

gsl_vector* getNodeValue()
{
    gsl_vector* nodeValue
        = gsl_vector_calloc((g_sliceNum + 1) * (g_sliceNum + 1));
    int isBdryNode = 0, i = 0, m = 0, n = 0;
    double temp = 0.0, h = 0.0;
    for (i = 0; i < (g_sliceNum + 1) * (g_sliceNum + 1); ++i) {
        m = i % (g_sliceNum + 1);
        n = i / (g_sliceNum + 1);
        isBdryNode = m == 0 || m == g_sliceNum || n == 0 || n == g_sliceNum;
        if (isBdryNode)
            gsl_vector_set(nodeValue, i, g_bdry((double)m * h, (double)n * h));
        else {
            temp = gsl_vector_get(
                g_innerNodeValue, (n - 1) * (g_sliceNum - 1) + m - 1);
            gsl_vector_set(nodeValue, i, temp);
        }
    }
    return nodeValue;
}

double getError(func u)
{
    double error = 0.0, h = 1.0 / (double)g_sliceNum;
    int i = 0, m = 0, n = 0;
    for (i = 0; i < (g_sliceNum - 1) * (g_sliceNum - 1); ++i) {
        m = i % (g_sliceNum - 1) + 1;
        n = i / (g_sliceNum - 1) + 1;
        error += (u((double)m * h, (double)n * h)
                     - gsl_vector_get(g_innerNodeValue, i))
            * (u((double)m * h, (double)n * h)
                  - gsl_vector_get(g_innerNodeValue, i));
    }
    return h * sqrt(error);
}

void sfemTest(void)
{
    int i = 0, j = 0;
    for (i = 0; i < (g_sliceNum - 1) * (g_sliceNum - 1); ++i) {
        for (j = 0; j < (g_sliceNum - 1) * (g_sliceNum - 1); ++j)
            printf("%f\t", gsl_spmatrix_get(g_stiffnessMatrix, i, j));
        printf("|\t%f\t|", gsl_vector_get(g_innerNodeValue, i));
        printf("|\t%f\n", gsl_vector_get(g_loadVector, i));
    }
}

void sfemFinal(void)
{
    gsl_spmatrix_free(g_stiffnessMatrix);
    gsl_spmatrix_free(g_stiffnessMatrixExt);
    gsl_vector_free(g_innerNodeValue);
    gsl_vector_free(g_loadVector);
}
