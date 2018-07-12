/* eval_int_normal.c ********************************************************************/

#include "eval_int_normal.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


/**
 * @title Evaluate Integrand for the case of a multivariate normal distribution
 * @param int n: size of the sample (i.e. number of rows of U)
 * @param int q: dimension of the problem (i.e. number of columns of U + 1)
 * @param double U: "matrix" (actually vector in C) consisting of uniforms (e.g. a Sobol pointset)
 * @param double a: vector of length q (lower limits)
 * @param double b: vector of length q (upper limits)
 * @param double C: "matrix" (actually vector in C), lower triangular Cholesky factor
 * @param double ONE: largest number x<1 such that x != 1
 * @param double ZERO: smallest number x>0 such that x!=0
 * @return double mean, mean of g(U) using antithetic variates
 *
 * @note a and b can have -Inf / +Inf entries. While pnorm(..) would give the correct result, it safes time to check each time if the argument is +/- Inf and setting the value to 1/0 rather than calling pnorm().
 *
 * @author Erik Hintz
 */
double eval_int_normal(int n, int q, double *U, double *a, double *b, double *C, double ONE, double ZERO)
{
    double y[q-1], d, f, cone, diff;
	/* the following variables (ending in -a) are used to store the corresponding antithetic values */
    double ya[q-1], da, fa, conea, diffa;
    
    double tmp;
	double mean = 0;
	int i, j, l;

	/* now that .Call is being used, U and C are vectors. To access their elements, we use the rule
	   U[s,k] = U[k*numrows+s] */

	/* note that C starts initializing at 0 */

	/* for each row of U */
	for(j=0; j < n; j++){
        
		/* initialize d,f */
        /* check if entry of a is -Inf (see note above) */
        if(a[0] == R_NegInf){
            d  = 0;
            da = 0;
        } else {
            d  = pnorm( a[0] / C[0], 0, 1, 1, 0);
            da = pnorm( a[0] / C[0], 0, 1, 1, 0);
        }
        /* check if entry of b is +Inf (see note above) */
        if(b[0]==R_PosInf){
            diff  = 1 - d;
            diffa = 1 - da;
        } else {
            diff  = pnorm( b[0] / C[0], 0, 1, 1, 0) - d;
            diffa = pnorm( b[0] / C[0], 0, 1, 1, 0) - da;
            
        }
        
        f = diff;
        fa = diffa;
        
		/* and then go through the columns */
		for(i=0; i< q-1; i++){

			/* U[i*n+j] corresponds to U[j,i] in the orginal matrix */
			tmp = d + U[i*n+j] * diff;

			/* check if too close to 1 or 0 */
			if(tmp > ONE){
				tmp = ONE;
			}
			if(tmp < ZERO){
				tmp = ZERO;
			}
			y[i] = qnorm(tmp, 0, 1, 1, 0);

			/* now the same for the antithetic value */
			tmp = da+ (1-U[i*n+j]) * diffa;

			if(tmp > ONE){
				tmp = ONE;
			}
			if(tmp < ZERO){
				tmp = ZERO;
			}
            ya[i] = qnorm(tmp, 0, 1, 1, 0);
            
        

			/* calculate the scalar product sum C[i,j]y[j] for j=1 to (i-1) */
			cone  = 0;
			conea = 0;

            for(l=0; l<(i+1); l++){
				/* C[l*q+i+1] corresponds to C[i+1,l] in the original matrix */
				cone  += y[l]  * C[l*q+i+1];
				conea += ya[l] * C[l*q+i+1];
			}
            
            /* calculate new d and diff = e-d. */
            if(a[i+1] == R_NegInf){
                d  = 0;
                da = 0;
            } else {
                d  = pnorm( (a[i+1] - cone) / C[(i+1)*(q+1)], 0, 1, 1, 0);
                da = pnorm( (a[i+1] - conea)/ C[(i+1)*(q+1)], 0, 1, 1, 0);
            }
            if(b[i+1] == R_PosInf){
                diff  = 1 - d;
                diffa = 1 - da;
            } else {
                diff  = pnorm( (b[i+1] - cone) /C[(i+1)*(q+1)], 0, 1, 1, 0) - d;
                diffa = pnorm( (b[i+1] - conea)/C[(i+1)*(q+1)], 0, 1, 1, 0) - da;
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
 * @title R Interface for eval_int_normal
 * @param ...same parameters as in eval_int_normal above
 * @return mean( f(U) ) where f is the integrand and U is the point-set
 * @author Erik Hintz
 */
SEXP eval_int_normal_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP ONE, SEXP ZERO)
{
	double res = eval_int_normal(INTEGER(n)[0], INTEGER(q)[0], REAL(U), REAL(a), REAL(b), REAL(C),REAL(ONE)[0], REAL(ZERO)[0]);
	return ScalarReal(res);
}
