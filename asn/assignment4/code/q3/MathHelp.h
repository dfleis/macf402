/*
	Header file for MathHelp
	Contains some useful functions for our purposes
*/
#ifndef MATHHELP_H
#define MATHHELP_H

#include <vector>

class MathHelp {
	public:
		static double stdNormPdf(double x);
		static double stdNormCdf(double x);
		static double mean(std::vector<double> v);
		static double sd(std::vector<double> v);
};
#endif