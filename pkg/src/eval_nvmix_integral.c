#include "eval_nvmix_integral.h"


/**
 * @title Evaluate Integrand for a Normal Variance Mixture
 * @param lower q-vector of lower evaluation limits
 * @param upper q-vector of upper evaluation limits
 * @param U vector representing '(n, q+1)-matrix' of uniforms
 *        (e.g. Sobol pointset) and realizations of sqrt(mix).
 *        The '(n, q+1)-matrix' is of the form [M1, U, M2], where the first
 *        column (M1) consists of realizations of sqrt(mix), the last column
 *        (M2) consists of the antithetic realizations of (M1) and the middle
 *        part is an (n, q-1)-matrix of uniforms (e.g. a Sobol point-set).
 * @param n sample size (i.e., number of rows of U)
 * @param q dimension minus 1 (i.e., number of columns of U minus 1)
 * @param cholScale lower triangular Cholesky factor as vector
 * @param ZERO smallest number x > 0 such that x != 0
 * @param ONE   largest number x < 1 such that x != 1
 * @return mean estimate mean(f(U)) of E(f(U)) using antithetic variates
 * @note lower and upper can have -/+Inf entries. While pnorm() would give the
         correct result, it safes time to check each time if the argument is -/+Inf
 *       and setting the value to 0/1 rather than calling pnorm().
 * @author Erik Hintz and Marius Hofert
 */
double eval_nvmix_integral_c(double *lower, double *upper, double *U, int n, int q,
                             double *cholScale, double ZERO, double ONE)
{
    double y[q-1],  sqrtmix,  d,  diff,  f,  scprod;
    double ya[q-1], sqrtmixa, da, diffa, fa, scproda; /* corresponding antithetic values */
    /* y:       vector to save phi^{-1}(dj+uj(ej-dj)) */
    /* sqrtmix: used to store sqrt(F_w^{-1}(u_0)) */
    /* d:       current values of di from the paper */
    /* diff:    current values of (ei-di) from the paper */
    /* f:       current value of (e1-d1) * (e2-d2) * ... * (ei-di) */
    /* scprod:  scalar product sum cholScale_{ij} y_j */
    double tmp; /* to store temporary values */
    double mean = 0; /* to store the result */
    int i, j, l; /* counters for loops */

    /* For each row of U (for access, use U[s,k] = U[k * numrows + s]) */
    for(j = 0; j < n; j++){

        /* Grab the realizations of sqrt(mix) = sqrt(W) */
        sqrtmix  = U[j];
        sqrtmixa = U[q*n+j];

        /* Check if entry of lower is -Inf */
        if(lower[0] == R_NegInf){
            d  = 0;
            da = 0;
        } else {
            d  = pnorm(lower[0] / (cholScale[0] * sqrtmix),  0, 1, 1, 0);
            da = pnorm(lower[0] / (cholScale[0] * sqrtmixa), 0, 1, 1, 0);
        }

        /* Check if entry of b is +Inf */
        if(upper[0] == R_PosInf){
            diff  = 1 - d;
            diffa = 1 - da;
        } else {
            diff  = pnorm(upper[0] / (cholScale[0] * sqrtmix),  0, 1, 1, 0) - d;
            diffa = pnorm(upper[0] / (cholScale[0] * sqrtmixa), 0, 1, 1, 0) - da;
        }

        /* Go through all q-1 columns (without first and last) */
        /* For better readability, we start at i = 0 */
        f  = diff;
        fa = diffa;
        for(i = 0; i < q-1; i++){
            /* U[i * n + j] corresponds to U[j,i] in the orginal matrix */
            tmp = d + U[(i+1) * n + j] * diff;

            /* Check if too close to 0 or 1 */
            if(tmp < ZERO){
                tmp = ZERO;
            }
            if(tmp > ONE){
                tmp = ONE;
            }
            y[i] = qnorm(tmp, 0, 1, 1, 0);

            /* The same for the antithetic value */
            tmp = da + (1-U[(i+1) * n + j]) * diffa;
            if(tmp < ZERO){
                tmp = ZERO;
            }
            if(tmp > ONE){
                tmp = ONE;
            }
            ya[i] = qnorm(tmp, 0, 1, 1, 0);

            /* Calculate the scalar product sum cholScale[i,j] y[j] for j = 1:(i-1) */
            scprod  = 0;
            scproda = 0;
            for(l = 0; l < (i+1); l++){
                /* cholScale[l * q + i+1] corresponds to cholScale[i+1,l] in the original matrix */
                scprod  += y[l]  * cholScale[l * q + i+1];
                scproda += ya[l] * cholScale[l * q + i+1];
            }

            /* Calculate new d and diff = e-d */
            if(lower[i+1] == R_NegInf){
                d  = 0;
                da = 0;
            } else {
                d  = pnorm((lower[i+1] / sqrtmix  - scprod)  / cholScale[(i+1)*(q+1)], 0, 1, 1, 0);
                da = pnorm((lower[i+1] / sqrtmixa - scproda) / cholScale[(i+1)*(q+1)], 0, 1, 1, 0);
            }
            if(upper[i+1] == R_PosInf){
                diff  = 1 - d;
                diffa = 1 - da;
            } else {
                diff  = pnorm((upper[i+1] / sqrtmix  - scprod)  / cholScale[(i+1)*(q+1)], 0, 1, 1, 0) - d;
                diffa = pnorm((upper[i+1] / sqrtmixa - scproda) / cholScale[(i+1)*(q+1)], 0, 1, 1, 0) - da;
            }
            f  *= diff;
            fa *= diffa;
        }
        mean += (f+fa)/2;
    }
    mean = mean/ n;
    return(mean);
}

/**
 * @title R Interface for eval_nvmix_integral_c()
 * @param see eval_nvmix_integral_c()
 * @return mean(f(U)) where f is the integrand and U specifies the point-set
 * @author Erik Hintz, Marius Hofert (polishing)
 */
SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP U, SEXP n, SEXP q,
			 SEXP cholScale, SEXP ZERO, SEXP ONE)
{
    double res = eval_nvmix_integral_c(REAL(lower), REAL(upper), REAL(U),
				       INTEGER(n)[0], INTEGER(q)[0],
				       REAL(cholScale), REAL(ZERO)[0],
                                       REAL(ONE)[0]);
    return ScalarReal(res);
}
