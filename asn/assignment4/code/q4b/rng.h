/*
	Header file for the random number generator
*/
#ifndef RNG_H
#define RNG_H

#include <vector>

class rng {
	public:
		static double runif();
		static double rnorm();
		static std::vector<double> rnorm(int n);
};
#endif