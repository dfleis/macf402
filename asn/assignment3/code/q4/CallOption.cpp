#include "CallOption.h"
#include <cmath>
#include "MathHelp.h"

CallOption::CallOption() {}
CallOption::CallOption(Asset A_in, double strike_in, double tau_in)
	:Option(A_in, strike_in, tau_in){} // call the parent constructor
double CallOption::payoff() {
	return fmax(A.getSpot() - strike, 0);
}
double CallOption::payoff(double some_spot) {
	return fmax(some_spot - strike, 0);
}
double CallOption::blackScholesPrice() {
	double R = A.getR();
	double spot = A.getSpot();
	
	return MathHelp::stdNormCdf(d1) * spot - MathHelp::stdNormCdf(d2) * strike * 
		exp(-R * tau);
}
double CallOption::blackScholesDelta() {
	return MathHelp::stdNormCdf(d1);
}
double CallOption::blackScholesDK() { // derivative wrt strike
	double r = A.getR();
	return -MathHelp::stdNormCdf(d2) * exp(-r * tau);
}
double CallOption::blackScholesDT() { // derivative wrt T
	double r = A.getR();
	double sigma = A.getSigma();
	
	return strike * exp(-r * tau) * (r * MathHelp::stdNormCdf(d2) + 
		sigma * MathHelp::stdNormPdf(d2) / (2 * sqrt(tau)) );
}