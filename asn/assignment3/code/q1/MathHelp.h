/*
	Header file for MathHelp
	Contains some useful functions for our purposes
*/
#ifndef MATHHELP_H
#define MATHHELP_H

#include <vector>

using namespace std;

class MathHelp {
	public:
		static double stdNormPdf(double x);
		static double stdNormCdf(double x);
		static double kernelSmoother(vector<double> h, vector<double> data_point);
		static double nadWatEstimator(std::vector<double> h_vec, 
			vector<double> x_unobs, vector<vector<double> > x_obs, vector<double> y_obs);
};
#endif