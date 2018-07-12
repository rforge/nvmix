//
//  eval_int_t.h
//

#ifndef eval_int_t_h
#define eval_int_t_h

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double eval_int_t(int n, int q, double *U, double *a, double *b,
             double *C, double nu, double ONE, double ZERO);

SEXP eval_int_t_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP nu, SEXP ONE, SEXP ZERO);

#endif /* eval_int_t_h */
