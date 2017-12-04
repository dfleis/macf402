#include "rng.h"
#include <stdlib.h>     /* srand, rand */
#include <cmath>

// constants for rnorm
const double a0 = 2.50662823884;
const double a1 = -18.61500062529;
const double a2 = 41.39119773534;
const double a3 = -25.44106049637;
const double a[] = {a0,a1,a2,a3};
const double b0 = -8.47351093090;
const double b1 = 23.08336743743;
const double b2 = -21.06224101826;
const double b3 = 3.13082909833;
const double b[] = {b0,b1,b2,b3};
const double c0 = 0.3374754822726147;
const double c1 = 0.9761690190917186;
const double c2 = 0.1607979714918209;
const double c3 = 0.0276438810333863;
const double c4 = 0.0038405729373609;
const double c5 = 0.0003951896511919;
const double c6 = 0.0000321767881768;
const double c7 = 0.0000002888167364;
const double c8 = 0.0000003960315187;
const double c[] = {c0,c1,c2,c3,c4,c5,c6,c7,c8};

double rng::runif() {
	return ((double) rand()) / RAND_MAX;
}

double rng::rnorm() {
	double v = runif();
	double u; 
	double z = 0;
	double A = 0; 
	double B = 0;
	
	u = v;
	if (v < 0.5) {
		u = 1 - v;
	} 
	
	if (u <= 0.92 & u >= 0.5) {
		for (int i = 0; i <= 3; i++) {
			A += a[i] * pow(u - 0.5, 2 * i + 1);
			B += b[i] * pow(u - 0.5, 2 * i + 2);
		}
		z = A / (1 + B);
	}	
	else {
		for (int i = 0; i <= 8; i++) {
			z += c[i] * pow(log( -log(1 - u) ), i);
		}
	}
	
	if (v < 0.5) return -z;
	else return z;
}
std::vector<double> rng::rnorm(int n) {
	std::vector<double> z(n);
	for (int i = 0; i < n; i++) {
		z.at(i) = rnorm();
	}
	return z;
}










































