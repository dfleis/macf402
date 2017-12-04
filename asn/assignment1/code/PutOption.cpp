#include "PutOption.h"
#include <cmath>

PutOption::PutOption() {}
PutOption::PutOption(Asset A_in, double strike_in, double tau_in)
	:Option(A_in, strike_in, tau_in){} // call the parent constructor
double PutOption::payoff() {
	return fmax(strike - A.getSpot(), 0);
}
double PutOption::payoff(double some_spot) {
	return fmax(strike - some_spot, 0);
}