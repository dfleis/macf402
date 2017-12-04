#include "MathHelp.h"
#include <cmath>

double MathHelp::log_stir(int n) { // use log(stirling's approx) to get a value of n! for large n
	if (n == 0) return 0; // we want Choose(N,0) = Choose(N,0) = 1
		// to do this requires 0! = 1, but we pass to log_fact the param. 0 in these border cases.
		// we know that log(0) is undefined, but we want exp(log(x)) = 1
		// so lets pretend log(0) = 0 since we later exp. it as exp(0) = 1
	return log(sqrt(2 * M_PI * n)) + n * (log(n) - 1);
}
double MathHelp::bin_coef(int n, int k) { 
	if (k > n) return 0; // this should never be called in our case

	if (n <= 1000) { // use this implementation for modest n in the bin. coef. Choose(n,k)
 		double r = 1; // result
 		
 		// since bin. coef. is symmetric, consider the half with small k for fewer iterations
 		if (k > n - k) k = n - k;
 		
		for (int i = 1; i <= k; i++) {
			r *= (n - i + 1)/((double)i);
		}
		return r;
	}
	else { // use this implementation for large n in the bin. coef. Choose(n,k);
    	   // apply log properties to log(n!/(k!(n-k)!)
		return (log_stir(n) - log_stir(k) - log_stir(n - k));
	}
}