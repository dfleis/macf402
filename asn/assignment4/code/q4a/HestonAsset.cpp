#include "HestonAsset.h"
#include <cmath> // pow

HestonAsset::HestonAsset() {}
HestonAsset::HestonAsset(double spot, double r, double vol, double sigma, double kappa, 
	double theta, double rho):Asset(spot, r) {
	
	this->vol = vol;
	this->sigma = sigma;
	this->kappa = kappa;
	this->theta = theta;
	this->rho = rho;
}
void HestonAsset::setVol(double vol) {
	this->vol = vol;
}
double HestonAsset::getVol() {
	return vol;
}
void HestonAsset::setSigma(double sigma) {
	this->sigma = sigma;
}
double HestonAsset::getSigma() {
	return sigma;
}
void HestonAsset::setKappa(double kappa) {
	this->kappa = kappa;
}
double HestonAsset::getKappa() {
	return kappa;
}
void HestonAsset::setTheta(double theta) {
	this->theta = theta;
}
double HestonAsset::getTheta() {
	return theta;
}
void HestonAsset::setRho(double rho) {
	this->rho = rho;
}
double HestonAsset::getRho() {
	return rho;
}
double HestonAsset::generateVolPath(double T, double z) {
	return hestonVol(vol, T, 1, z);
}
std::vector<double> HestonAsset::generateVolPath(double T, int n, std::vector<double> z) {
	// to do, check if sizes match
	std::vector<double> v(n + 1);
	v.at(0) = vol;
	
	for (int i = 1; i < v.size(); i++) {
		v.at(i) = hestonVol(v.at(i - 1), T, n, z.at(i - 1));
	}
	return v;
}
double HestonAsset::generatePricePath(double T, double z1, double z2) {
	return hestonPrice(spot, vol, T, 1, z1, z2);
}
std::vector<double> HestonAsset::generatePricePath(double T, int n, 
	std::vector<double> z1, std::vector<double> z2) { // to do, check if size of z matches
	std::vector<double> S(n + 1), v(n + 1);
	S.at(0) = spot;
	v = generateVolPath(T, n, z1);
	
	for (int i = 1; i < S.size(); i++) {
		S.at(i) = hestonPrice(S.at(i - 1), v.at(i - 1), T, n, z1.at(i - 1), z2.at(i - 1));
	}
	return S;
}
double HestonAsset::hestonVol(double v, double T, int n, double z) {
	double v_out = v + kappa * (theta - v) * (T/n) + sigma * sqrt(v * (T/n)) * z;
	if (v_out < 0) return 0;
	else return v_out;
}
double HestonAsset::hestonPrice(double S, double v, double T, int n, double z1, double z2) {
	double x = log(S);
	
	return exp(x + (r - 0.5 * v) * (T/n) + sqrt(v) * (rho * sqrt(T/n) * z1 + 
		sqrt(1 - pow(rho, 2)) * sqrt(T/n) * z2));
}




































