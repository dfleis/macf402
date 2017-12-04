/*
	Header file for an Option
	An option is an object that contains an asset, a strike, and a time to maturity
	Maybe it makes sense to implement this as an abstract class? 
	It doesn't make sense to be able to instantiate a general 'Option'
*/
#include "Asset.h"

#ifndef OPTION_H
#define OPTION_H
class Option {
	protected:
		Asset A; // underlying asset
		double strike; // strike price
		double tau;  // time to expiry
				
	public:
		Option(); // default constructor
		Option(Asset A_in, double strike_in, double tau_in);
		void setAsset(Asset A_in);
		Asset getAsset();
		void setStrike(double strike_in);
		double getStrike();
		void setTau(double tau_in);
		double getTau();
		
		virtual double payoff() = 0;
		virtual double payoff(double some_strike) = 0;
};
#endif