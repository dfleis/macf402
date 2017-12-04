#include "AriAsianOption.h"
#include <cmath>

AriAsianOption::AriAsianOption() {}
AriAsianOption::AriAsianOption(Asset A, double K, double t, double T, int n, char type)
	:AsianOption(A, K, t, T, n, type){} // call the parent constructor
double AriAsianOption::payoff() {
	double avg;
	double df = exp(-A.getR() * (T - t));
	double sum = 0;

	for (int i = 1; i < path.size(); i++) {
		sum += path.at(i);
	}
	avg = sum/n;

	if (type == 'c') return df * fmax(avg - K, 0);
	else if (type == 'p') return df * fmax(K - avg, 0);
	else return nan("1");
}
double AriAsianOption::price() {
	return nan("1"); // maybe numeric implementation w/ default control variate algo
}