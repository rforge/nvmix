//
//  evalfbonly.h
//  
//
//  Created by Erik Hintz on 6/7/18.
//
//

#ifndef evalfbonly_h
#define evalfbonly_h

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double evalfbonly(int n, int q, double *U, double *b,
                  double *C, double nu, double ONE, double ZERO);
SEXP evalfbonly_(SEXP n, SEXP q, SEXP U, SEXP b, SEXP C, SEXP nu, SEXP ONE, SEXP ZERO);


#endif /* evalfbonly_h */
