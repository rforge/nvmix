//
//  evalf.h
//  
//
//  Created by Erik Hintz on 6/7/18.
//
//

#ifndef evalf_h
#define evalf_h

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double evalf(int n, int q, double *U, double *a, double *b,
             double *C, double nu, double ONE, double ZERO);

SEXP evalf_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP nu, SEXP ONE, SEXP ZERO);

#endif /* evalf_h */
