#include "GeoAsianOption.h"
#include <cmath>
#include "MathHelp.h"

GeoAsianOption::GeoAsianOption() {}
GeoAsianOption::GeoAsianOption(Asset A, double K, double t, double T, int n, char type)
	:AsianOption(A, K, t, T, n, type){} // call the parent constructor
double GeoAsianOption::payoff() {
	double g_avg;
	double df = exp(-A.getR() * (T - t));
	double prod = 1;

	for (int i = 1; i < path.size(); i++) {
		prod *= path.at(i);
	}
	g_avg = pow(prod, 1.0/n);

	if (type == 'c') return df * fmax(g_avg - K, 0);
	else if (type == 'p') return df * fmax(K - g_avg, 0);
	else return nan("1");
}
double GeoAsianOption::price() {
	double S = A.getSpot();
	double r = A.getR();
	double mu_hat = computeMuHat();
	double d1 = computeD1();
	double d2 = computeD2();
	
	if (type == 'c') {
		return exp(-r * (T - t)) * (S * exp(mu_hat * (T - t)) * MathHelp::stdNormCdf(d1) - 
			K * MathHelp::stdNormCdf(d2));
	} else if (type == 'p') {
		return 0;
	}
	else return nan("1");
}
double GeoAsianOption::computeSigmaHat() {
	double sigma = A.getSigma();
	
	return sqrt( pow(sigma, 2) * (n + 1) * (2 * n + 1) / (6 * pow(n, 2)) );
}
double GeoAsianOption::computeMuHat() {
	double r = A.getR();
	double sigma = A.getSigma();
	double sigma_hat = computeSigmaHat();
	
	return (r - 0.5 * pow(sigma, 2)) * (n + 1) / (2 * n) + 0.5 * pow(sigma_hat, 2);
}
double GeoAsianOption::computeD1() {
	double S = A.getSpot();
	double r = A.getR();
	double sigma = A.getSigma();
	
	double sigma_hat = computeSigmaHat();
	double mu_hat = computeMuHat();
	
	return (log(S / K) + (mu_hat + 0.5 * pow(sigma_hat, 2)) * (T - t)) / (sigma_hat * sqrt(T - t));
}
double GeoAsianOption::computeD2() {
	double sigma = A.getSigma();
	double sigma_hat = computeSigmaHat();
	
	return computeD1() - sigma_hat * sqrt(T - t);
}






