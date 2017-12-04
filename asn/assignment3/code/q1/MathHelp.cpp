#include "MathHelp.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

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
double MathHelp::kernelSmoother(vector<double> h, vector<double> data_point) {
	// normal density kernel smoother for arbitrary dimensional point
	int dim = data_point.size();
	
	double K = 1/(h.at(0) * sqrt(2 * M_PI)) * exp(-1/2.0 * 
		pow(data_point.at(0)/h.at(0), 2));
	
	for (int i = 1; i < dim; i++) {
		K *= 1/(h.at(i) * sqrt(2 * M_PI)) * exp(-1/2.0 * pow(data_point.at(i)/h.at(i), 2));
	}

	return K;
}
double MathHelp::nadWatEstimator(vector<double> h_vec, vector<double> x_unobs, 
	vector<vector<double> > x_obs, vector<double> y_obs) {
	
	if (h_vec.size() != x_unobs.size()) {
		cout << "WARNING in nad_wat_estimator: Dimensions to not agree" << endl;
	}
	int n = x_unobs.size(); // number of x variables to smooth over
		
	double num = 0;
	double den = 0;
	
	vector<double> data_point(n); // stores n dimensional difference of unobs - obs
		
	for (int i = 0; i < x_obs.at(0).size(); i++) { // loop over observed data
		for (int j = 0; j < n; j++) { // loop over dimensions to store difference
			data_point.at(j) = x_unobs.at(j) - x_obs.at(j).at(i);
		}	

		num += kernelSmoother(h_vec, data_point) * y_obs.at(i);
		den += kernelSmoother(h_vec, data_point);
	}
	
	return num/den;
}




































