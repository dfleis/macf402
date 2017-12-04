/*
	Header file for an Asset
	An asset is an object with a spot, interest rate, and (optional) volatility 
*/
#ifndef ASSET_H
#define ASSET_H

#include <vector>

class Asset {
	private:
		double spot; // asset spot value
		double r; // cont. comp. interest rate
		double sigma; // volatility parameter
		
		double lognormalModel(double S, double T, int n, double Z);
		
	public:
		Asset(); // default constructor
		Asset(double spot, double r, double sigma);
		void setSpot(double spot);
		double getSpot();
		void setR(double r);
		double getR();
		void setSigma(double sigma);
		double getSigma();
		
		std::vector<double> generatePath(double T, int n, std::vector<double> z);
};
#endif
