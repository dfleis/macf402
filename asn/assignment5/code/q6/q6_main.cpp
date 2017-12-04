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
double bs_price_call(double S, double K, double r, double sigma, double dt) {
	if (dt == 0) 
		return fmax(S - K, 0);
	double d1 = computeD1(S, K, r, sigma, dt);
	double d2 = computeD2(S, K, r, sigma, dt);
	return stdNormCdf(d1) * S - stdNormCdf(d2) * K * exp(-r * dt);
}

int main() {
	srand(402);
	int n_min = 2; // log10 min to compute
	int n_max = 9; // log10 max to compute
	
	double S0 = 100;
	double K = 100;
	double mu = 0.02;
	double sigma = 0.5;
	double T = 1.0;
	double h = 5.0/365;
	double z_alpha = 1.645;
	string filepath = "../../data/q6_dat.csv";
	
	double C0 = bs_price_call(S0, K, mu, sigma, T);
	double d1 = computeD1(S0, K, mu, sigma, T);
	double delta0 = stdNormCdf(d1) * S0;
	double var = z_alpha * sigma * sqrt(h) * delta0;
	
	// initialize & open csv
	ofstream outfile;
	outfile.open(filepath);
	// write header
	outfile << "N,short,short_se, long, long_se" << endl;
	
	for (int n = n_min; n <= n_max; n++) {
		int N = pow(10, n);
		vector<double> z = rng::rnorm(N);
		double St, Ct, dC;
	
		int short_sum = 0;
		int long_sum = 0;

		// compute var breaks
		for (int i = 0; i < N; i++) { // loop over each asset simulation
			St = S0 * exp((mu - 0.5 * pow(sigma, 2)) * h + sigma * sqrt(h) * z.at(i));
			Ct = bs_price_call(St, K, mu, sigma, T - h);
			dC = Ct - C0;
			if (dC > var) 
				short_sum++;
			if (dC < -var)
				long_sum++;
		}
		
		// compute proportions
		double short_pct = short_sum / (double)N;
		double long_pct = long_sum / (double)N;
		
		// compute SE
		double short_se = sqrt( short_pct * (1 - short_pct) / N );
		double long_se = sqrt( short_pct * (1 - short_pct) / N );

		// write data to csv
		outfile << N << "," << short_pct << "," << short_se << "," << 
			long_pct << "," << long_se << endl;
	}
	outfile.close();
}









































