#include "Option.h"
#include <iostream>
#include <cmath>
#include "MathHelp.h"

Option::Option() {}
Option::Option(Asset A_in, double strike_in, double tau_in) {
	A = A_in;
	strike = strike_in;
	tau = tau_in;
	
	double spot = A.getSpot();
	double sigma = A.getSigma();
	double R = A.getR();
	
	d1 = computeD1(spot, strike, tau, R, sigma);
	d2 = computeD2(spot, strike, tau, R, sigma);
}	
void Option::updateOption() {
	double spot = A.getSpot();
	double sigma = A.getSigma();
	double R = A.getR();
		
	d1 = computeD1(spot, strike, tau, R, sigma);
	d2 = computeD2(spot, strike, tau, R, sigma);	
}	
void Option::setAsset(Asset A_in) {
	A = A_in;
}
Asset Option::getAsset() {
	return A;
}
void Option::setStrike(double strike_in) {
	strike = strike_in;
}
double Option::getStrike() {
	return strike;
}
void Option::setTau(double tau_in) {
	tau = tau_in; 
}
double Option::getTau() {
	return tau;
}
void Option::setD1(double d1_in) {
	d1 = d1_in; 
}
double Option::getD1() {
	return d1;
}
void Option::setD2(double d2_in) {
	d2 = d2_in; 
}
double Option::getD2() {
	return d2;
}
double Option::computeD1(double S, double K, double tau, double R, double sigma) {
	return 1/(sigma * sqrt(tau)) * (log(S/K) + (R + 0.5 * pow(sigma,2)) * tau);
}
double Option::computeD2(double S, double K, double tau, double R, double sigma) {
	return computeD1(S, K, tau, R, sigma) - sigma * sqrt(tau);
}

double Option::blackScholesVega() {
	if (A.getSigma() == 0) {
		std::cout << "WARNING in black_scholes_vega(): sigma = 0, is this correct?" << std::endl;
		return 0;
	}
	double spot = A.getSpot();	
	return spot * MathHelp::stdNormPdf(d1) * sqrt(tau);
}
double Option::blackScholesD2K2() { // second derivative wrt strike
	double sigma = A.getSigma();
	double r = A.getR();
	
	return MathHelp::stdNormPdf(d2) * exp(-r * tau) / ( strike * sigma * sqrt(tau));
}

double Option::impliedVolNR(double V_obs, double x_0, double eps) {
	if (V_obs != V_obs) { // exploit the fact that nan is not equal to itself
		return std::nan("1"); // this is for simplicity when dealing with options w/o asks
	}
	
	double sigma_guess = x_0;
	
	setAsset(Asset(A.getSpot(), A.getR(), sigma_guess));
	updateOption();
	
	double V_guess = blackScholesPrice();
			
	while(fabs(V_guess - V_obs) >= eps) {
		double f = blackScholesPrice() - V_obs;
		double f_prime = blackScholesVega();
		
		if (f_prime == 0) {
			std::cout << "WARNING in newton_raphson: computed vega = 0, divide by zero" << std::endl;
			break;
		}
				
		sigma_guess = sigma_guess - f/f_prime;
		
		setAsset(Asset(A.getSpot(), A.getR(), sigma_guess));
		updateOption();
		
		V_guess = blackScholesPrice();
	}
	return sigma_guess;	
}
double Option::impliedVolBisect(double V_obs, double a, double b, double eps) {
	if (V_obs != V_obs) { // exploit the fact that nan is not equal to itself
		return std::nan("1"); // this is for simplicity when dealing with options w/o asks
	}

	double sigma_guess = (a + b)/2;

	setAsset(Asset(A.getSpot(), A.getR(), sigma_guess));
	updateOption();
	
	double V_guess = blackScholesPrice();
		
	while(fabs(V_guess - V_obs) >= eps) {
		if (V_guess > V_obs) { // guessed price too high, vol too high
			b = sigma_guess; // shrink right endpoint
			sigma_guess = (a + b)/2;
		}
		else { // guessed price too low, vol too low
			a = sigma_guess; // shrink left endpoint
			sigma_guess = (a + b)/2;
		}

		setAsset(Asset(A.getSpot(), A.getR(), sigma_guess));
		updateOption();
		
		V_guess = blackScholesPrice();
		
	}	
	return sigma_guess;
}




































