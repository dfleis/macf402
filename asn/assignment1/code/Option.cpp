#include "Option.h"

Option::Option() {}
Option::Option(Asset A_in, double strike_in, double tau_in) {
	A = A_in;
	strike = strike_in;
	tau = tau_in;
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