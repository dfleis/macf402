/*
	Header file for a VanillaOption
	A vanilla option is a child of an option that contains a specific implementation
	of a payoff function
*/
#ifndef VANILLAOPTION_H
#define VANILLAOPTION_H

#include "Option.h"

class VanillaOption: public Option {
	private:
		double computeD1();
		double computeD2();
		
	public:
 		VanillaOption(); // default constructor
 		VanillaOption(Asset A, double K, double t, double T, char type);
		
		double payoff();
		double blackScholesPrice();
		double blackScholesDelta();
};
#endif