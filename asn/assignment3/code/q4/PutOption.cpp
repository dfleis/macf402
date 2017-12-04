#include "PutOption.h"
#include <cmath>
#include "MathHelp.h"

PutOption::PutOption() {}
PutOption::PutOption(Asset A_in, double strike_in, double tau_in)
	:Option(A_in, strike_in, tau_in){} // call the parent constructor
double PutOption::payoff() {
	return fmax(strike - A.getSpot(), 0);
}
double PutOption::payoff(double some_spot) {
	return fmax(strike - some_spot, 0);
}
double PutOption::blackScholesPrice() {
	double R = A.getR();
	double spot = A.getSpot();
	
	return -MathHelp::stdNormCdf(-d1) * spot + MathHelp::stdNormCdf(-d2) * strike * 
		exp(-R * tau);
}
double PutOption::blackScholesDelta() {
	return -MathHelp::stdNormCdf(-d1);
}
double PutOption::blackScholesDK() { // derivative wrt strike
	double r = A.getR();
	return MathHelp::stdNormCdf(-d2) * exp(-r * tau);
}
double PutOption::blackScholesDT() { // derivative wrt T
	double r = A.getR();
	double sigma = A.getSigma();
	
	return strike * exp(-r * tau) * (-r * MathHelp::stdNormCdf(-d2) + 
		sigma * MathHelp::stdNormPdf(d2) / (2 * sqrt(tau)) );
}