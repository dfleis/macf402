#include "AsianOption.h"

AsianOption::AsianOption() {}
AsianOption::AsianOption(Asset A, double K, double t, double T, int n, char type)
	:Option(A, K, t, T, type) {
	this->n = n;
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