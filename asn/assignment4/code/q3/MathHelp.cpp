#include "MathHelp.h"
#include <cmath>

double MathHelp::stdNormPdf(double x) {
	return 1/sqrt(2 * M_PI) * exp(-1/2.0 * pow(x, 2));
}
double MathHelp::stdNormCdf(double x_in) {
	// we notice that the approximation is much better for positive x
  	// use symmetry of the normal CDF to compute the negative x values
  	double x = fabs(x_in);
  
	double b0 = 0.2316419;
  	double b1 = 0.31938530;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double t = 1/(1 + b0 * x);
  
	double y = 1 - stdNormPdf(x) * (b1*t + b2*pow(t,2) + b3*pow(t,3) + b4*pow(t,4) + b5*pow(t,5));
	
	if (x_in < 0) {
		y = 1 - y;
	}
	return y;
}
double MathHelp::mean(std::vector<double> v) {
	double sum = 0;
	for (int i = 0; i < v.size(); i++) {
		sum += v.at(i);
	}
	return sum / v.size();
}
double MathHelp::sd(std::vector<double> v) {
	double xbar = mean(v);
	double sq_sum = 0;
	for (int i = 0; i < v.size(); i++) {
		sq_sum += pow((v.at(i) - xbar), 2);
	}
	return sqrt( sq_sum / (v.size() - 1) );
}