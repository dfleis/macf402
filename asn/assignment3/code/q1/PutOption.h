/*
	Header file for a PutOption
	A put option is a child of an option that contains a specific implementation
	of a payoff function
*/
#ifndef PUTOPTION_H
#define PUTOPTION_H

#include "Option.h"

class PutOption: public Option {
	public:
 		PutOption(); // default constructor
 		PutOption(Asset A_in, double strike_in, double tau_in);
		
		double payoff();
		double payoff(double some_spot);
		double blackScholesPrice();
		double blackScholesDelta();
};
#endif