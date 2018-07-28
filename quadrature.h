#ifndef QUADRATURE_H
#define QUADRATURE_H

typedef double (*FuncWithArg)(const double, const double, const void*);

double getQuadrature(const FuncWithArg f, const void *arg);

#endif