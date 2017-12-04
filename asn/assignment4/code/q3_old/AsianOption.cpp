#include "AsianOption.h"
#include <cmath> // fmax
#include <stdlib.h> // exit

AsianOption::AsianOption() {}
AsianOption::AsianOption(Asset A, double K, double t, double T, int n, char type) {
	this->A = A;
	this->K = K;
	this->t = t;
	this->T = T;
	this->n = n;
	this->type = type;
}
Asset AsianOption::getAsset() {
	return A;
}
void AsianOption::setAsset(Asset A) {
	this->A = A;
}
double AsianOption::getTStart() {
	return t;
}
void AsianOption::setTStart(double t) {
	this->t = t;
}
double AsianOption::getTEnd() {
	return T;
}
void AsianOption::setTEnd(double T) {
	this->T = T;
}
double AsianOption::getNMonitoringPoints() {
	return n;
}
void AsianOption::setNMonitoringPoints(int n) {
	this->n = n;
}
std::vector<double> AsianOption::getPath() {
	return path;
}
void AsianOption::setPath(std::vector<double> path) {
	this->path = path;
}