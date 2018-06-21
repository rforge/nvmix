/* evalf.c ********************************************************************/

#include "evalf.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


/**
 * @title Evaluate Integrand If Not All Lower Endpoints Are -Inf
 * @param <TODO: describe all parameters (each on its own line)>
 * @return <TODO>
 * @author Erik Hintz
 */
double evalf(int n, int q, double *U, double *a, double *b,
             double *C, double nu, double ONE, double ZERO)
{
	double finit, dinit;
	double y[q-1], sumysq, d, f, cone, ctwo, diff, tmp;
	/* the following variables (ending in -a) are used to store the corresponding antithetic values */
	double ya[q-1], sumysqa, da, fa, conea, ctwoa, diffa, tmpa;
	double mean=0;
	int i, j, l;

	/* now that .Call is being used, U and C are vectors. To access their elements, we use the rule
	   U[s,k] = U[k*numrows+s] */

	/* note that C starts initializing at 0 */

	/* initialize d and f. These values can be reused */
	dinit = pt( a[0] / C[0], nu, 1, 0);
	finit = pt( b[0] / C[0], nu, 1, 0) - dinit;

	/* for each row of U */
	for(j=0; j < n; j++){

		/* initialize y */

		/* U[j] corresponds to U[j,0] in the orginal matrix */
		y[0]  = qt(dinit + U[j]*finit, nu, 1, 0);
		ya[0] = qt(dinit + (1-U[j])*finit, nu, 1, 0);

		sumysq  = y[0] * y[0];
		sumysqa = ya[0] * ya[0];

		/* C[1] corresponds to C[1,0] in the original matrix */
		cone  = y[0] * C[1];
		conea = ya[0] * C[1];

		/* C[q+1] corresponds to C[1,1] in the original matrix */
		ctwo  = sqrt( (nu+1)/(nu+sumysq) ) / C[q+1];
		ctwoa = sqrt( (nu+1)/(nu+sumysqa) ) / C[q+1];

		d  = pt( ( a[1] - cone) * ctwo, nu+1, 1, 0);
		da = pt( ( a[1] - conea) * ctwoa, nu+1, 1, 0);

		diff  = pt( ( b[1] - cone) * ctwo, nu+1, 1, 0) - d;
		diffa = pt( ( b[1] - conea) * ctwoa, nu+1, 1, 0) - da;

		f  = diff * finit;
		fa = diffa * finit;

		/* and then go through the columns */
		/* first column has already been taken care of */

		for(i=1; i< (q-1); i++){

			/* U[i*n+j] corresponds to U[j,i] in the orginal matrix */
			tmp = d+ U[i*n+j]*diff;

			/* check if too close to 1 or 0 */
			if(tmp > ONE){
				tmp = ONE;
			}
			if(tmp < ZERO){
				tmp = ZERO;
			}

			y[i] = qt(tmp, nu+i, 1, 0) * sqrt( (nu + sumysq) / (nu+i));
			sumysq += y[i] * y[i];


			/* now the same for the antithetic value */
			tmpa = da+ (1-U[i*n+j])*diffa;

			if(tmpa > ONE){
				tmpa = ONE;
			}
			if(tmpa < ZERO){
				tmpa = ZERO;
			}

			ya[i] = qt(tmpa, nu+i, 1, 0) * sqrt( (nu + sumysqa) / (nu+i));
			sumysqa += ya[i] * ya[i];

			/* calculate the scalar product */
			cone  = 0;
			conea = 0;

			for(l=0; l<(i+1); l++){
				/* C[l*q+i+1] corresponds to C[i+1,l] in the original matrix */
				cone  += y[l]* C[l*q+i+1];
				conea += ya[l]* C[l*q+i+1];
			}

			/* C[(i+1)*(q+1)] corresponds to C[i+1,i+1] in the original matrix */
			ctwo  = sqrt( (nu+i+1) / (nu + sumysq)) / C[(i+1)*(q+1)];
			ctwoa = sqrt( (nu+i+1) / (nu + sumysqa)) / C[(i+1)*(q+1)];

			d  = pt( (a[i+1] - cone)*ctwo, nu+i+1, 1, 0);
			da = pt( (a[i+1] - conea)*ctwoa, nu+i+1, 1, 0);
			diff  = pt( (b[i+1] - cone)*ctwo, nu+i+1, 1, 0) - d;
			diffa = pt( (b[i+1] - conea)*ctwoa, nu+i+1, 1, 0) - da;

			f  *= diff;
			fa *= diffa;
		}
		mean += (f+fa)/2;
	}
	mean = mean/ n;
	return(mean);
}

/**
 * @title <TODO: provide title here>
 * @param <TODO: describe all parameters (each on its own line)>
 * @return <TODO>
 * @author Erik Hintz
 */
SEXP evalf_(SEXP n, SEXP q, SEXP U, SEXP a, SEXP b, SEXP C, SEXP nu, SEXP ONE, SEXP ZERO)
{
	double res = evalf(INTEGER(n)[0], INTEGER(q)[0], REAL(U), REAL(a), REAL(b), REAL(C),
			   REAL(nu)[0], REAL(ONE)[0], REAL(ZERO)[0]);
	return ScalarReal(res);
}
