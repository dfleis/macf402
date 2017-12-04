#include "CallOption.h"
#include <cmath>

CallOption::CallOption() {}
CallOption::CallOption(Asset A_in, double strike_in, double tau_in)
	:Option(A_in, strike_in, tau_in){} // call the parent constructor
double CallOption::payoff() {
	return fmax(A.getSpot() - strike, 0);
}
double CallOption::payoff(double some_spot) {
	return fmax(some_spot - strike, 0);
}
