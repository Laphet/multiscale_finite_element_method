#include "spfem.h"

static int g_sliceNum;
static coefficient g_A;
static func g_f, g_f1, g_f2;
static gsl_spmatrix* g_stiffnessMatrix;
static gsl_spmatrix* g_stiffnessMatrixExt;
static gsl_vector* g_periodicNodeValue;
static gsl_vector* g_loadVector;

static int g_m, g_n, g_flag;
static const int NUM_OF_INTE = 10, NUM_OF_LOAD = 4;
static const double ERR_TOL = 1.0e-7;
static const int MAX_ITER = 1024;

void sfemInit(int sliceNum)
{
    g_sliceNum = sliceNum;
    g_stiffnessMatrix
        = gsl_spmatrix_alloc(g_sliceNum * g_sliceNum, g_sliceNum * g_sliceNum);
    g_stiffnessMatrixExt = gsl_spmatrix_alloc(
        g_sliceNum * g_sliceNum, (g_sliceNum + 1) * (g_sliceNum + 1));
    g_periodicNodeValue = gsl_vector_calloc(g_sliceNum * g_sliceNum);
    g_loadVector = gsl_vector_calloc(g_sliceNum * g_sliceNum);
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

static double getLocalStiffness(m, n, flag)
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

    for (i = 0; i < g_sliceNum * g_sliceNum; ++i) {
        m = i % g_sliceNum;
        n = i / g_sliceNum;
        m_ = (m - 1) < 0 ? g_sliceNum : m - 1;
        n_ = (n - 1) < 0 ? g_sliceNum : n - 1;
        temp = getLocalStiffness(m, n, 1) + getLocalStiffness(m, n_, 2)
            + getLocalStiffness(m_, n_, 3) + getLocalStiffness(m, n_, 4);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m, n), temp);
        temp = getLocalStiffness(m, n_, 7) + getLocalStiffness(m, n, 5);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n), temp);
        temp = getLocalStiffness(m, n, 8) + getLocalStiffness(m_, n, 6);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m, n + 1), temp);
        temp = getLocalStiffness(m_, n, 5) + getLocalStiffness(m_, n_, 7);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m_, n), temp);
        temp = getLocalStiffness(m_, n_, 6) + getLocalStiffness(m, n_, 8);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m, n_), temp);
        temp = getLocalStiffness(m, n, 0);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n + 1), temp);
        temp = getLocalStiffness(m_, n, 9);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m_, n + 1), temp);
        temp = getLocalStiffness(m_, n_, 0);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m_, n_), temp);
        temp = getLocalStiffness(m, n_, 9);
        gsl_spmatrix_set(
            g_stiffnessMatrixExt, i, getIndexInVector(m + 1, n_), temp);
    }
    for (i = 0; i < g_sliceNum * g_sliceNum; ++i) {
        m = i % g_sliceNum;
        n = i / g_sliceNum;
        for (j = 0; j < i + 1; ++j) {
            m_ = j % g_sliceNum;
            n_ = j / g_sliceNum;
            temp = gsl_spmatrix_get(
                g_stiffnessMatrixExt, i, getIndexInVector(m_, n_));
            gsl_spmatrix_set(g_stiffnessMatrix, i, j, temp);
            gsl_spmatrix_set(g_stiffnessMatrix, j, i, temp);
        }
    }
}

static void setLoadVector(int isWithDivF)
{
    int m = 0, n = 0, m_ = 0, n_ = 0, i = 0;
    double temp = 0.0;

    for (i = 0; i < g_sliceNum * g_sliceNum; ++i) {
        m = i % g_sliceNum;
        n = i / g_sliceNum;
        m_ = (m - 1 < 0) ? g_sliceNum : m - 1;
        n_ = (n - 1 < 0) ? g_sliceNum : n - 1;
        temp = getLocalLoad(m, n, 0, isWithDivF)
            + getLocalLoad(m_, n, 1, isWithDivF)
            + getLocalLoad(m_, n_, 2, isWithDivF)
            + getLocalLoad(m, n_, 3, isWithDivF);
        gsl_vector_set(g_loadVector, i, temp);
    }
}

int solvePeriodicPDE(coefficient A, func f)
{
    g_A = A;
    g_f = f;
    setStiffnessMatrix();
    setLoadVector(0);
    const gsl_splinalg_itersolve_type* T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve* work = gsl_splinalg_itersolve_alloc(
        T, (g_sliceNum - 1) * (g_sliceNum - 1), 0);
    int status = 0, iter = 0;
    do {
        status = gsl_splinalg_itersolve_iterate(g_stiffnessMatrix,
            g_loadVector, ERR_TOL, g_periodicNodeValue, work);
        ++iter;
    } while (status == GSL_CONTINUE && iter < MAX_ITER);
    gsl_splinalg_itersolve_free(work);
    return 1;
}

int solvePeriodicPDEwithDivF(coefficient A, func f1, func f2)
{
    g_A = A;
    g_f1 = f1;
    g_f2 = f2;
    setLoadVector(1);
    const gsl_splinalg_itersolve_type* T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve* work = gsl_splinalg_itersolve_alloc(
        T, (g_sliceNum - 1) * (g_sliceNum - 1), 0);
    int status = 0, iter = 0;
    do {
        status = gsl_splinalg_itersolve_iterate(g_stiffnessMatrix,
            g_loadVector, ERR_TOL, g_periodicNodeValue, work);
        ++iter;
    } while (status == GSL_CONTINUE && iter < MAX_ITER);
    gsl_splinalg_itersolve_free(work);
    return 1;
}

gsl_vector* getNodeValue()
{
    gsl_vector* nodeValue
        = gsl_vector_calloc((g_sliceNum + 1) * (g_sliceNum + 1));
    int i = 0, m = 0, n = 0;
    double temp = 0.0;
    for (i = 0; i < (g_sliceNum + 1) * (g_sliceNum + 1); ++i) {
        m = i % (g_sliceNum + 1);
        n = i / (g_sliceNum + 1);
        m = (m == g_sliceNum) ? 0 : m;
        n = (n == g_sliceNum) ? 0 : n;
        temp = gsl_vector_get(g_periodicNodeValue, n * g_sliceNum + m);
        gsl_vector_set(nodeValue, i, temp);
    }

    return nodeValue;
}

double getError(func u)
{
    double error = 0.0, h = 1.0 / (double)g_sliceNum;
    int i = 0, m = 0, n = 0;
    for (i = 0; i < g_sliceNum * g_sliceNum; ++i) {
        m = i % g_sliceNum;
        n = i / g_sliceNum;
        error += (u((double)m * h, (double)n * h)
                     - gsl_vector_get(g_periodicNodeValue, i))
            * (u((double)m * h, (double)n * h)
                  - gsl_vector_get(g_periodicNodeValue, i));
    }
    return h * sqrt(error);
}

void sfemTest(void)
{
    int i = 0, j = 0;
    for (i = 0; i < g_sliceNum * g_sliceNum; ++i) {
        for (j = 0; j < g_sliceNum * g_sliceNum; ++j)
            printf("%f\t", gsl_spmatrix_get(g_stiffnessMatrix, i, j));
        printf("|\t%f\t|", gsl_vector_get(g_periodicNodeValue, i));
        printf("|\t%f\n", gsl_vector_get(g_loadVector, i));
    }
}

void sfemFinal(void)
{
    gsl_spmatrix_free(g_stiffnessMatrix);
    gsl_spmatrix_free(g_stiffnessMatrixExt);
    gsl_vector_free(g_periodicNodeValue);
    gsl_vector_free(g_loadVector);
}
