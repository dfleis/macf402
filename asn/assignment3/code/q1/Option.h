/*
	Header file for an Option
	An option is an object that contains an asset, a strike, and a time to maturity
	Maybe it makes sense to implement this as an abstract class? 
	It doesn't make sense to be able to instantiate a general 'Option'
*/
#ifndef OPTION_H
#define OPTION_H

#include "Asset.h"
#include <cmath>

class Option {
	protected:
		Asset A; // underlying asset
		double strike; // strike price
		double tau;  // time to expiry
		double d1;
		double d2;
				
	public:
		Option(); // default constructor
		Option(Asset A_in, double strike_in, double tau_in);
		void updateOption();
		void setAsset(Asset A_in);
		Asset getAsset();
		void setStrike(double strike_in);
		double getStrike();
		void setTau(double tau_in);
		double getTau();
		void setD1(double d1_in);
		double getD1();
		void setD2(double d2_in);
		double getD2();
				
		double computeD1(double S, double K, double tau, double R, double sigma);
		double computeD2(double S, double K, double tau, double R, double sigma);
		double blackScholesVega();
		double impliedVolNR(double V_obs, double x_0, double eps = pow(10,-6));
		double impliedVolBisect(double V_obs, double a, double b, double eps = pow(10,-6));
						
		virtual double payoff() = 0;
		virtual double payoff(double some_strike) = 0;
		virtual double blackScholesPrice() = 0;
		virtual double blackScholesDelta() = 0;
};
#endif