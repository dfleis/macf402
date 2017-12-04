#include "VanillaOption.h"
#include <cmath> // fmax
#include <stdlib.h> // exit
#include "MathHelp.h"

VanillaOption::VanillaOption() {}
VanillaOption::VanillaOption(Asset A, double K, double t, double T, char type)
	:Option(A, K, t, T, type){} // call the parent constructor
double VanillaOption::payoff() {
	if (type == 'c') return fmax(A.getSpot() - K, 0);
	else if (type == 'p') return fmax(K - A.getSpot(), 0);
	else exit(EXIT_FAILURE);
}
double VanillaOption::blackScholesPrice() {
	double R = A.getR();
	double S = A.getSpot();
	double d1 = computeD1();
	double d2 = computeD2();
	
	if (type == 'c') {
		return MathHelp::stdNormCdf(d1) * S - MathHelp::stdNormCdf(d2) * K * 
			exp(-R * (T - t));
	} else if (type == 'p') {
		return -MathHelp::stdNormCdf(-d1) * S + MathHelp::stdNormCdf(-d2) * K * 
		exp(-R * (T - t));
	} else exit(EXIT_FAILURE);
}
double VanillaOption::blackScholesDelta() {
	double d1 = computeD1();

	if (type == 'c') return MathHelp::stdNormCdf(d1);
	else if (type == 'p') return -MathHelp::stdNormCdf(-d1);
	else exit(EXIT_FAILURE);
}
double VanillaOption::computeD1() {
	double S = A.getSpot();
	double R = A.getR();
	double sigma = A.getSigma();
	
	return 1/(sigma * sqrt(T - t)) * (log(S/K) + (R + 0.5 * pow(sigma,2)) * (T - t));
}
double VanillaOption::computeD2() {
	double sigma = A.getSigma();
	return computeD1() - sigma * sqrt(T - t);
}