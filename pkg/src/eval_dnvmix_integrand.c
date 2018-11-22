#include "eval_dnvmix_integrand.h"


/**
 * @title Evaluate Integrand for a the log-density of a normal variance mixture
 * @param W current_n-vector of realizations of W (sorted)
 * @param maha2_2 n-vector of mahalabonis distances (maha^2/2) (sorted)
 * @param current_n RQMC sample size to estimate the density (length of W)
 * @param n number of evaluation points (length of maha2_2)
 * @param d dimension of the nvmix distribution
 * @param lrdet log(sqrt(det(scale)); 'scale' is the scale matrix of the dist'n
 * @return ldensities n-vector of estimated log-densities
 * @note ...
 * @author Erik Hintz a
 */
double eval_dnvmix_integrand_c(double *W, double *maha2_2, int current_n, int n,
                               int d, double lrdet)
{
    double ldensities[n]; /* result vector */
    double c[n]; /* vector to store c_i */
    double current_maha, current_c, next_c, current_W, c_max, sum_expc;
    int maxindex, found_max;
    int i, j, l; /* counters for loops */
    int startindex = 0; /* C starts counting at 0 */
    
    /* Some constants that we can re-use: */
    neglogcurrent_n = - log(current_n);
    d2 = -d/2;
    
    /* For each evaluation point in maha2_2 */
    for(j = 0; j < n; j++){
        /* Grab current maha2_2 */
        current_maha = maha2_2[j];
        
        /* Calculate c_i starting from i = startindex until max reached */
        i = startindex;
        found_max = 0;
        /* first c_i: */
        current_W = W[i];
        current_c = d2 * log(2*pi*current_W) - lrdet - current_maha/current_W;
        /* save this so that we can re-use it later. To save memory, we always
         fill the vector c from the beginning, no matter what. */
        c[i - startindex] = current_c;
        i += 1;
        while(found_max = 0){
            current_W = W[i];
            next_c = d2 * log(2*pi*current_W) - lrdet - current_maha/current_W;
            c[i - startindex] = next_c;
            /* Did we find the maximum? c is first increasing, then decreasing
             as a function of W */
            if(next_c < current_c){
                found_max = 1;
                c_max = current_c;
                /* The maximum occured in the previous index */
                maxindex = i - 1;
            }
            current_c = next_c;
            i += 1;
        }
        
        /* Calculate sum_1^current_n exp(ci-cmax):
        /* Term for i = max is exp(0) = 1 */
        sum_expc = 1;
        for(l = 0, l < current_n, l++){
            /* Was c_l already calculated? */
            if( l >= startindex && l <= (maxindex + 1) && !(l==maxindex)){
                sum_expc += exp(c[l - startindex] - c_max);
            } else {
                current_W = W[l];
                sum_expc += exp( d2 * log(2*pi*current_W) - lrdet - current_maha/current_W - cmax);
            }
        }
        
        /* the position of c_max is increasing in maha2_2. Since the latter is
         sorted in increasing order, we use startindex = maxindex for the
         next maha2_2 value */
        startindex = maxindex;
        
        /* Done: */
        ldensities[j] = neglogcurrent_n + c_max + log(sum_expc);
    }
    
    return(ldensities);
}


/**
 * @title R Interface for eval_dnvmix_integrand_c()
 * @param see eval_dnvmix_integrand_c() above
 * @return see eval_dnvmix_integrand_c() above
 * @author Erik Hintz
 */


SEXP eval_dnvmix_integrand(SEXP W, SEXP maha2_2, SEXP current_n, SEXP n, SEXP d,
			 SEXP lrdet)
{
    double res = eval_dnvmix_integrand_c(REAL(W), REAL(maha2_2),
				       INTEGER(current_n)[0], INTEGER(n)[0], INTEGER(d)[0],
				       REAL(lrdet));
    return Real(res);
}
