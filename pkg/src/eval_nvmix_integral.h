#ifndef eval_nvmix_integral_h
#define eval_nvmix_integral_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


double eval_nvmix_integral_c(double *a, double *b, double *U, int n, int q,
			     double *C, double ZERO, double ONE);
SEXP eval_nvmix_integral(SEXP a, SEXP b, SEXP U, SEXP n, SEXP q, SEXP C,
			 SEXP ZERO, SEXP ONE);

#endif /* eval_nvmix_integral_h */
