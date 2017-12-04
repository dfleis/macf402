#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "rng.h"

using namespace std;

double stdNormPdf(double x) {
	return 1/sqrt(2 * M_PI) * exp(-1/2.0 * pow(x, 2));
}
double stdNormCdf(double x_in) {
  	double x = fabs(x_in);
  
	double b0 = 0.2316419;
  	double b1 = 0.31938530;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double t = 1/(1 + b0 * x);
  
	double y = 1 - stdNormPdf(x) * (b1*t + b2*pow(t,2) + b3*pow(t,3) + b4*pow(t,4) + b5*pow(t,5));
	
	if (x_in < 0) {
		y = 1 - y;
	}
	return y;
}
double computeD1(double S, double K, double r, double sigma, double dt) {
	return (log(S / K) + (r + 0.5 * pow(sigma, 2)) * dt) / (sigma * sqrt(dt));
}
double computeD2(double S, double K, double r, double sigma, double dt) {
	return computeD1(S, K, r, sigma, dt) - sigma * sqrt(dt);
}
double bs_price_put(double S, double K, double r, double sigma, double dt) {
	if (dt == 0) 
		return fmax(S - K, 0);
	double d1 = computeD1(S, K, r, sigma, dt);
	double d2 = computeD2(S, K, r, sigma, dt);
	return stdNormCdf(-d2) * K * exp(-r * dt) - stdNormCdf(-d1) * S;
}

int main() {
	srand(402);
	
	// parameters
	double S0 = 100;
	double K = 100;
	double r = 0.03;
	double sigma = 0.2;
	double t_max = 1;
	double S_max = 3 * K;
	
	int N = 10;
	int J = 5;
	
	// increments
	double dS = S_max / J;
	double dt = t_max / N;
	double nu = dt / pow(dS, 2);
	
	// create the tridiagonal matrix values
	vector<double> aa(J - 1), bb(J - 1), cc(J - 1);
	for (int j = 1; j < J; j++) {
		aa.at(j - 1) =  (-0.5 * r * dt * j + 0.5 * pow(sigma, 2) * dt * pow(j, 2)) / (1 + r * dt);
		bb.at(j - 1) = (1 - pow(sigma, 2) * dt * pow(j, 2)) / (1 + r * dt);
		cc.at(j - 1) = (0.5 * r * dt * j + 0.5 * pow(sigma, 2) * dt * pow(j, 2)) / (1 + r * dt);
	}
	
	// initialize matrix of solutions
	vector<vector<double> > V;
	// matrix.resize(num_of col, vector<double> (num_of_row , init_value));
	V.resize(N + 1, vector<double> (J + 1, 0));
	
	// place initial & boundary conditions
	for (int j = 0; j < V.at(0).size(); j++) { // initial conditions
		V.at(N).at(j) = fmax(K - j * dS, 0);
	}
	for (int n = 0; n < V.size(); n++) { // boundary conditions
		V.at(n).at(0) = exp(-r * (t_max - n * dt)) * K;
	}
	
	// FTCS algorithm
	for (int n = V.size() - 2; n > -1; n--) {
		for (int j = 1; j < V.at(0).size() - 1; j++) {
			V.at(n).at(j) = V.at(n + 1).at(j - 1) * aa.at(j - 1) + 
				V.at(n + 1).at(j) * bb.at(j - 1) + 
				V.at(n + 1).at(j + 1) * cc.at(j - 1);
		}
	}
	
	// write data
	string filepath = "../../data/q5_dat_N" + to_string(N) + "_J" + 
		to_string(J) + ".csv";	
	
	ofstream outfile;
	outfile.open(filepath);
	
	for (int i = 0; i < V.size(); i++) {
		for (int j = 0; j < V.at(0).size() - 1; j++) {
			outfile << V.at(i).at(j) << ",";
		}
		outfile << V.at(i).at(V.at(0).size() - 1) << endl;
	}
	outfile.close();
}










































