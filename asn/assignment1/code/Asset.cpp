#include "Asset.h"

Asset::Asset() {}
Asset::Asset(double spot_in, double R_in) {
	spot = spot_in;
	R = R_in;
}
Asset::Asset(double spot_in, double R_in, double sigma_in) {
	spot = spot_in;
	R = R_in;
	sigma = sigma_in;
}
void Asset::setSpot(double spot_in) {
	spot = spot_in;
}
double Asset::getSpot() {
	return spot;
}
void Asset::setR(double R_in) {
	R = R_in;
}
double Asset::getR() {
	return R;
}
void Asset::setSigma(double sigma_in) {
	sigma = sigma_in;
}
double Asset::getSigma() {
	return sigma;
}