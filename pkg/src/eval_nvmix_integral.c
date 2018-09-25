#include "eval_nvmix_integral.h"


/**
 * @title Evaluate Integrand for a Normal Variance Mixture
 * @param a vector of length q (lower limits)
 * @param b vector of length q (upper limits)
 * @param U vector representing '(n, q+1)-matrix' of uniforms
 *        (e.g. Sobol pointset) and realizations of sqrt(mix).
 *        The '(n, q+1)-matrix' is of the form [M1, U, M2], where the first
 *        column (M1) consists of realizations of sqrt(mix), the last column
 *        (M2) consists of the antithetic realizations of (M1) and the middle
 *        part is an (n, q-1)-matrix of uniforms (e.g. a Sobol point-set).
 * @param n sample size (i.e., number of rows of U)
 * @param q dimension minus 1 (i.e., number of columns of U minus 1)
 * @param C lower triangular Cholesky factor as vector
 * @param ZERO smallest number x > 0 such that x != 0
 * @param ONE   largest number x < 1 such that x != 1
 * @return mean estimate of E(g(U)) using antithetic variates
 * @note a and b can have -/+Inf entries. While pnorm() would give the correct
 *       result, it safes time to check each time if the argument is -/+Inf
 *       and setting the value to 0/1 rather than calling pnorm().
 * @author Erik Hintz, Marius Hofert (polishing)
 */
double eval_nvmix_integral_c(double *a, double *b, double *U, int n, int q,
			     double *C, double ZERO, double ONE)
{
    double y [q-1], sqrtmix,  d,  f,  cone,  diff; /* TODO provide brief description for each variable */
    double ya[q-1], sqrtmixa, da, fa, conea, diffa; /* corresponding antithetic values */
    double tmp; /* TODO provide brief description for each variable */
    double mean = 0; /* TODO provide brief description for each variable */
    int i, j, l; /* TODO provide brief description for each variable */

    /* For each row of U (for access, use U[s,k] = U[k * numrows + s]) */
    for(j=0; j < n; j++){

        /* Grab the realizations of sqrt(mix) = sqrt(W) TODO: correct ('sqrt(W)')? */
        sqrtmix  = U[j];
        sqrtmixa = U[q*n+j];

        /* Check if entry of a is -Inf */
        if(a[0] == R_NegInf){
            d  = 0;
            da = 0;
        } else {
            d  = pnorm(a[0] / (C[0] * sqrtmix),  0, 1, 1, 0);
            da = pnorm(a[0] / (C[0] * sqrtmixa), 0, 1, 1, 0);
        }

        /* Check if entry of b is +Inf */
        if(b[0] == R_PosInf){
            diff  = 1 - d;
            diffa = 1 - da;
        } else {
            diff  = pnorm(b[0] / (C[0] * sqrtmix),  0, 1, 1, 0) - d;
            diffa = pnorm(b[0] / (C[0] * sqrtmixa), 0, 1, 1, 0) - da;
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

	    /* Calculate the scalar product sum C[i,j] y[j] for j = 1:(i-1) */
	    cone  = 0;
	    conea = 0;
            for(l = 0; l < (i+1); l++){
		/* C[l * q + i+1] corresponds to C[i+1,l] in the original matrix */
		cone  += y[l]  * C[l * q + i+1];
		conea += ya[l] * C[l * q + i+1];
	    }

            /* Calculate new d and diff = e-d */
            if(a[i+1] == R_NegInf){
                d  = 0;
                da = 0;
            } else {
                d  = pnorm((a[i+1] / sqrtmix  - cone)  / C[(i+1)*(q+1)], 0, 1, 1, 0);
                da = pnorm((a[i+1] / sqrtmixa - conea) / C[(i+1)*(q+1)], 0, 1, 1, 0);
            }
            if(b[i+1] == R_PosInf){
                diff  = 1 - d;
                diffa = 1 - da;
            } else {
                diff  = pnorm((b[i+1] / sqrtmix  - cone)  / C[(i+1)*(q+1)], 0, 1, 1, 0) - d;
                diffa = pnorm((b[i+1] / sqrtmixa - conea) / C[(i+1)*(q+1)], 0, 1, 1, 0) - da;
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
SEXP eval_nvmix_integral(SEXP a, SEXP b, SEXP U, SEXP n, SEXP q, SEXP C,
			 SEXP ZERO, SEXP ONE)
{
    double res = eval_nvmix_integral_c(REAL(a), REAL(b), REAL(U), INTEGER(n)[0],
				       INTEGER(q)[0], REAL(C), REAL(ZERO)[0],
				       REAL(ONE)[0]);
    return ScalarReal(res);
}
