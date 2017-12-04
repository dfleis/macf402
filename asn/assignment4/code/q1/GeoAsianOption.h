/*
	Header file for a GeoAsianOption
*/
#ifndef GEOASIANOPTION_H
#define GEOASIANOPTION_H

#include "AsianOption.h"
#include <vector>

class GeoAsianOption: public AsianOption {
	private:		
		double computeSigmaHat();
		double computeMuHat();	
		double computeD1();
		double computeD2();
		
	public:
 		GeoAsianOption(); // default constructor
		GeoAsianOption(Asset A, double K, double t, double T, int n, char type);
		
		double payoff();	
		double price();
};
#endif