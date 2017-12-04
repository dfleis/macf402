#include "Option.h"

Option::Option() {}
Option::Option(Asset A, double K, double t, double T, char type) {
	this->A = A;
	this->K = K;
	this->t = t;
	this->T = T;
	this->type = type;
}		
void Option::setAsset(Asset A) {
	this->A = A;
}
Asset Option::getAsset() {
	return A;
}
void Option::setStrike(double K) {
	this->K = K;
}
double Option::getStrike() {
	return K;
}
void Option::setTStart(double t) {
	this->t = t;
}
double Option::getTStart() {
	return t;
}
void Option::setTEnd(double T) {
	this->T = T;
}
double Option::getTEnd() {
	return T;
}
void Option::setType(char type) {
	this->type = type;
}
char Option::getType() {
	return type;
}
































