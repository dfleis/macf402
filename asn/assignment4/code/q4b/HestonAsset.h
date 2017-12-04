/*
	Header file for an HestonAsset
*/
#ifndef HESTONASSET_H
#define HESTONASSET_H

#include "Asset.h"
#include <vector>

class HestonAsset: public Asset {
	private:
		double vol; // vol 
		double sigma; // vol of vol
		double kappa; // reversion rate 
		double theta; // long run variance 
		double rho;
		
		double hestonVol(double v, double T, int n, double z);
		double hestonPrice(double S, double v, double T, int n, double z1, double z2);
		
	public:
		HestonAsset(); // default constructor
		HestonAsset(double spot, double r, double vol, double sigma, double kappa, 
			double theta, double rho); // heston
		void setVol(double vol);
		double getVol();
		void setSigma(double sigma);
		double getSigma();
		void setKappa(double kappa);
		double getKappa();
		void setTheta(double theta);
		double getTheta();
		void setRho(double rho);
		double getRho();
		
		double generateVolPath(double T, double z);
		std::vector<double> generateVolPath(double T, int n, std::vector<double> z);
		double generatePricePath(double T, double z1, double z2);
		std::vector<double> generatePricePath(double T, int n, 
			std::vector<double> z1, std::vector<double> z2);
};
#endif
