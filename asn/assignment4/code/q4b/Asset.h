/*
	Header file for an Asset
	An asset is an object with a spot, interest rate, and (optional) volatility 
*/
#ifndef ASSET_H
#define ASSET_H

#include <vector>

class Asset {
	protected:
		double spot; // asset spot value
		double r; // cont. comp. interest rate
		double sigma; // vol (lognormal)
		
		double lognormalModel(double S, double T, int n, double Z);
		
	public:
		Asset(); // default constructor
		virtual ~Asset(){};
		Asset(double spot, double r); // constructor w/o vol
		Asset(double spot, double r, double sigma); 
		void setSpot(double spot);
		double getSpot();
		void setR(double r);
		double getR();
		void setSigma(double sigma);
		double getSigma();
		
		double generatePath(double T, double z);
		std::vector<double> generatePath(double T, int n, std::vector<double> z);
};
#endif
