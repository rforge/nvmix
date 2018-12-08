#ifndef eval_dnvmix_integrand_h
#define eval_dnvmix_integrand_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void eval_dnvmix_integrand_c(double *W, double *maha2_2, int current_n, int n,
                               int d, int k, double lrdet, double *ldensities);
SEXP eval_dnvmix_integrand(SEXP W, SEXP maha2_2, SEXP current_n, SEXP n, SEXP d,
                           SEXP k, SEXP lrdet);

#endif /* eval_dnvmix_integrand_h */
