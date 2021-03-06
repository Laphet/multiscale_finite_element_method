#ifndef SFEM_H
#define SFEM_H

#include "gsl/gsl_splinalg.h"
#include "gsl/gsl_spmatrix.h"
#include "gsl/gsl_vector.h"
#include "numint.h"
#include <math.h>
#include <stdio.h>

typedef struct coefficient {
    func a1, a2, a3;
} coefficient;

void sfemInit(int sliceNum);

int solvePDE(coefficient A, func f, func bdry);

int solvePDEwithDivF(coefficient A, func f1, func f2, func bdry);

void sfemTest(void);

double getError(func u);

void sfemFinal(void);

#endif
