/*
	Header file for an Option
	An option is an object that contains an asset, a strike, and a time to maturity
	Maybe it makes sense to implement this as an abstract class? 
	It doesn't make sense to be able to instantiate a general 'Option'
*/
#ifndef OPTION_H
#define OPTION_H

#include "Asset.h"

class Option {
	protected:
		Asset A; // underlying asset
		double K; // strike price
		double t; // current time
		double T; // expiry time
		char type; // call or put
				
	public:
		Option(); // default constructor
		Option(Asset A, double K, double t, double T, char type);
		void setAsset(Asset A);
		Asset getAsset();
		void setStrike(double K);
		double getStrike();
		void setTStart(double t);
		double getTStart();
		void setTEnd(double T);
		double getTEnd();
		void setType(char type);
		char getType();
								
		virtual double payoff() = 0;
};
#endif