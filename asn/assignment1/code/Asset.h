/*
	Header file for an Asset
	An asset is an object with a spot, interest rate, and (optional) volatility 
*/
#ifndef ASSET_H
#define ASSET_H
class Asset {
	private:
		double spot; // asset spot value
		double R; // cont. comp. interest rate
		double sigma; // volatility parameter
		
	public:
		Asset(); // default constructor
		Asset(double spot_in, double R_in);
		Asset(double spot_in, double R_in, double sigma_in);
		void setSpot(double spot_in);
		double getSpot();
		void setR(double R_in);
		double getR();
		void setSigma(double sigma_in);
		double getSigma();
};
#endif