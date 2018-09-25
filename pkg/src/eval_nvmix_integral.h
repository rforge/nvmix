#ifndef eval_nvmix_integral_h
#define eval_nvmix_integral_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


double eval_nvmix_integral_c(int n, int q, double *U, double *a, double *b,
			     double *C, double ONE, double ZERO);
SEXP eval_nvmix_integral(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C,
			 SEXP ONE, SEXP ZERO);

#endif /* eval_nvmix_integral_h */
