#include "Asset.h"
#include <cmath> // pow

Asset::Asset() {}
Asset::Asset(double spot, double r, double sigma) {
	this->spot = spot;
	this->r = r;
	this->sigma = sigma;
}
void Asset::setSpot(double spot) {
	this->spot = spot;
}
double Asset::getSpot() {
	return spot;
}
void Asset::setR(double r) {
	this->r = r;
}
double Asset::getR() {
	return r;
}
void Asset::setSigma(double sigma) {
	this->sigma = sigma;
}
double Asset::getSigma() {
	return sigma;
}
std::vector<double> Asset::generatePath(double T, int n, std::vector<double> z) {
	// to do, check if size of z matches
	
	std::vector<double> S(n + 1);
	S.at(0) = spot;
	
	for (int i = 1; i < S.size(); i++) {
		S.at(i) = lognormalModel(S.at(i - 1), T, n, z.at(i - 1));
	}
	return S;
}
double Asset::lognormalModel(double S, double T, int n, double Z) {
	return S * exp( (r - 0.5 * pow(sigma, 2)) * T/n + sigma * sqrt(T/n) * Z );
}














