/*
	Header file for a AriAsianOption
*/
#ifndef ARIASIANOPTION_H
#define ARIASIANOPTION_H

#include "AsianOption.h"
#include <vector>

class AriAsianOption: public AsianOption {
	public:
 		AriAsianOption(); // default constructor
		AriAsianOption(Asset A, double K, double t, double T, int n, char type);
		
		double payoff();
		double price();
};
#endif