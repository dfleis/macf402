/*
	Header file for a CallOption
	A call option is a child of an option that contains a specific implementation
	of a payoff function
*/
#ifndef CALLOPTION_H
#define CALLOPTION_H

#include "Option.h"

class CallOption: public Option {
	public:
 		CallOption(); // default constructor
 		CallOption(Asset A_in, double strike_in, double tau_in);
		
		double payoff();
		double payoff(double some_spot);
		double blackScholesPrice();
		double blackScholesDelta();
};
#endif