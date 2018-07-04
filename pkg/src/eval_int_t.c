/* eval_int_t.c ********************************************************************/

#include "eval_int_t.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


/**
 * @title Evaluate Integrand for the case of a multivariate t distribution
 * @param int n: size of the sample (i.e. number of rows of U)
 * @param int q: dimension of the problem (i.e. number of columns of U)
 * @param double U: "matrix" (actually vector in C) consisting of uniforms (e.g. a Sobol pointset)
 * @param double a: vector of length q (lower limits)
 * @param double b: vector of length q (upper limits)
 * @param double C: "matrix" (actually vector in C), lower triangular Cholesky factor
 * @param double nu: degrees of freedom
 * @param double ONE: largest number x<1 such that x != 1
 * @param double ZERO: smallest number x>0 such that x!=0
 * @return double mean, mean of g(U) using antithetic variates
 *
 * @note a and b can have -Inf / +Inf entries. While pnorm(..) would give the correct result, it safes time to check each time if the argument is +/- Inf and setting the value to 1/0 rather than calling pnorm().
 *
 * @author Erik Hintz
 */
double eval_int_t(int n, int q, double *U, double *a, double *b,
             double *C, double nu, double ONE, double ZERO)
{
	double y[q-1], d, f, cone, diff, ssqrtnu;
	/* the following variables (ending in -a) are used to store the corresponding antithetic values */
    double ya[q-1], da, fa, conea, diffa, ssqrtnua;
    
    double tmp;
	double mean = 0;
	int i, j, l;

	/* now that .Call is being used, U and C are vectors. To access their elements, we use the rule
	   U[s,k] = U[k*numrows+s] */

	/* note that C starts initializing at 0 */

	/* for each row of U */
	for(j=0; j < n; j++){
        
        /* calculate X/sqrt(nu) where X is chi_(nu), simulated from the first dimension of U */
        ssqrtnu  = sqrt( qchisq( U[j], nu, 1, 0)   / nu);
        ssqrtnua = sqrt( qchisq( 1-U[j], nu, 1, 0) / nu);
        
		/* initialize d,f */
        /* check if entry if a is -Inf (see note above) */
        if(a[0] == R_NegInf){
            d  = 0;
            da = 0;
        } else {
            d  = pnorm( ssqrtnu  * a[0] / C[0], 0, 1, 1, 0);
            da = pnorm( ssqrtnua * a[0] / C[0], 0, 1, 1, 0);
        }
        /* check if entry if b is +Inf (see note above) */
        if(b[0]==R_PosInf){
            diff  = 1 - d;
            diffa = 1 - da;
        } else {
            diff  = pnorm( ssqrtnu  * b[0] / C[0], 0, 1, 1, 0) - d;
            diffa = pnorm( ssqrtnua * b[0] / C[0], 0, 1, 1, 0) - da;
            
        }
        
        f = diff;
        fa = diffa;
        
		/* and then go through the remaining columns */
        /* we start at i=0 to make it  more readable */
		for(i=0; i< q-1; i++){

			/* U[(i+1)*n+j] corresponds to U[j,i+1] in the orginal matrix */
			tmp = d+ U[(i+1)*n+j] * diff;

			/* check if too close to 1 or 0 */
			if(tmp > ONE){
				tmp = ONE;
			}
			if(tmp < ZERO){
				tmp = ZERO;
			}
			y[i] = qnorm(tmp, 0, 1, 1, 0);

			/* now the same for the antithetic value */
			tmp = da+ (1-U[(i+1)*n+j]) * diffa;

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
                d  = pnorm( (ssqrtnu  * a[i+1] - cone) / C[(i+1)*(q+1)], 0, 1, 1, 0);
                da = pnorm( (ssqrtnua * a[i+1] - conea)/ C[(i+1)*(q+1)], 0, 1, 1, 0);
            }
            if(b[i+1] == R_PosInf){
                diff  = 1 - d;
                diffa = 1 - da;
            } else {
                diff  = pnorm( (ssqrtnu  * b[i+1] - cone) /C[(i+1)*(q+1)], 0, 1, 1, 0) - d;
                diffa = pnorm( (ssqrtnua * b[i+1] - conea)/C[(i+1)*(q+1)], 0, 1, 1, 0) - da;
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
 * @title R Interface for eval_int_t
 * @param ...same parameters as in eval_int_t above
 * @return mean( f(U) ) where f is the integrand and U is the point-set
 * @author Erik Hintz
 */
SEXP eval_int_t_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP nu, SEXP ONE, SEXP ZERO)
{
	double res = eval_int_t(INTEGER(n)[0], INTEGER(q)[0], REAL(U), REAL(a), REAL(b), REAL(C),
			   REAL(nu)[0], REAL(ONE)[0], REAL(ZERO)[0]);
	return ScalarReal(res);
}
