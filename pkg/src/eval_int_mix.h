//
//  eval_int_mix.h
//

#ifndef eval_int_mix_h
#define eval_int_mix_h

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double eval_int_mix(int n, int q, double *U, double *a, double *b,
             double *C, double ONE, double ZERO);

SEXP eval_int_mix_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP ONE, SEXP ZERO);

#endif /* eval_int_mix_h */
