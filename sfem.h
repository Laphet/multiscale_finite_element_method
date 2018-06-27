#ifndef SFEM_H
#define SFEM_H

#include "gsl/gsl_spmatrix.h"
#include "gsl/gsl_splinalg.h"
#include "gsl/gsl_vector.h"
#include "numint.h"
#include <stdio.h>

typedef struct coefficient
{
    func a_1, a_2, a_3;
} coefficient;

double zero(double x, double y) { return 0.0; }

void sfemInit(int sliceNum);

void setCoefficient(coefficient A);

void solveEPDE(func f, func bdry);

void solveEPDEwithDivF(func f1, func f2, func bdry);

void sfemFinal(void);

#endif