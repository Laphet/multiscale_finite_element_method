#ifndef SPFEM_H
#define SPFEM_H

#include "gsl/gsl_splinalg.h"
#include "gsl/gsl_spmatrix.h"
#include "gsl/gsl_vector.h"
#include "numint.h"
#include <math.h>
#include <stdio.h>

typedef struct coefficient {
    func a1, a2, a3;
} coefficient;

void spfemInit(int sliceNum);

int solvePeriodicPDE(coefficient A, func f);

int solvePeriodicPDEwithDivF(coefficient A, func f1, func f2);

void spfemTest(void);

double getError(func u);

void spfemFinal(void);

#endif
