#ifndef SFEM_H
#define SFEM_H

#include "gsl/gsl_spmatrix.h"
#include "gsl/gsl_splinalg.h"
#include "gsl/gsl_vector.h"
#include "numint.h"
#include <stdio.h>

typedef struct coefficient
{
    func a1, a2, a3;
}coefficient;

double zero(double x, double y);

void sfemInit(int sliceNum);

void setCoefficient(coefficient A);

int solvePDE(func f, func bdry);

void solvePDEwithDivF(func f1, func f2, func bdry);

void sfemFinal(void);

#endif
