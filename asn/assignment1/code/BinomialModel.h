/*
	Header file for BinomialModel
	A binomial model is an operation (or set of operations) on an object and 
	thus is a static class that cannot be instantiated
*/
#include "Asset.h"
#include "Option.h"
#include "MathHelp.h"
#include <cmath>

#ifndef BINOMIALMODEL_H
#define BINOMIALMODEL_H
class BinomialModel {
	public:
		static double p; // risk neutral prob. of up state
		/* 
			Are factors really static attributes of BinomialModels?
			On the one hand it doesn't make sense to have BinomialModel be instantiated, 
			on the other hand up and down factors will vary across assets.
			Furthermore, up and down factors clearly should not be attributes of an Asset
		*/ 
		static double u; // up factor
		static double d; // down factor
		
		// price asset given n up steps given N total steps
		static double assetPrice(Asset A, double tau, int n, int N);
		// value option back to time = 0 given N step binomial model
		static double optionValuation(Option* O, int N);
		// back out implied vol. estimates via bisection algorithm
		static double impliedVol(double V_obs, Option* O, int N, double eps = pow(10,-6));
		
};
#endif