//
//  eval_int_normal.h
//

#ifndef eval_int_normal_h
#define eval_int_normal_h

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double eval_int_normal(int n, int q, double *U, double *a, double *b,
             double *C, double ONE, double ZERO);

SEXP eval_int_normal_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP ONE, SEXP ZERO);

#endif /* eval_int_normal_h */
