#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>
#include "rng.h"
#include "MathHelp.h"
#include "HestonAsset.h"
#include "VanillaOption.h"

using namespace std;

vector<double> CI(vector<double> v, double Z_a_2) {
	vector<double> conf_int(3);
	int n = v.size();
	double avg = MathHelp::mean(v);
	double stddev = MathHelp::sd(v);
		
	conf_int.at(0) = avg - Z_a_2 * stddev / sqrt(n);
	conf_int.at(1) = avg;
	conf_int.at(2) = avg + Z_a_2 * stddev / sqrt(n);
	
	return conf_int;
}
void writeMatrix(vector<vector<double> > data, string filepath, string header) {
	ofstream outfile;
	outfile.precision(6);
	
	outfile.open(filepath);
	outfile << header << endl;
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < data.at(0).size(); j++) {
			outfile << data.at(i).at(j) << ",";
		}
		outfile << endl;
	}
	outfile.close();
}
double conditionalXi(double rho, vector<double> v, double dt, vector<double> z) {
	double xi = 0;
	for (int i = 0; i < v.size(); i++) {
		xi += -0.5 * pow(rho, 2) * v.at(i) * dt + rho * sqrt(v.at(i) * dt) * z.at(i);
	}
	return exp(xi);
}
double conditionalSigma(vector<double> v, double dt) {
	double T = v.size() * dt;
	double sigma = 0;
	for (int i = 0; i < v.size(); i++) {
		sigma += v.at(i) * dt;
	}
	return sqrt(sigma/T);
}
double conditionalD1(double S, double K, double r, double sigma, double rho, double xi, double T) { 
	double num = log( S * xi / K ) + (r + 0.5 * pow(sigma, 2) * (1 - pow(rho, 2))) * T;
	double den = sigma * sqrt(T * (1 - pow(rho, 2)));
	return num/den;
}
double conditionalD2(double S, double K, double r, double sigma, double rho, double xi, double T) { 
	return conditionalD1(S, K, r, sigma, rho, xi, T) - sigma * sqrt(T * (1 - pow(rho, 2)));
}
vector<double> conditionalMC(HestonAsset HA, double K, double T, int n, int N) {
	double rho = HA.getRho();
	double S0 = HA.getSpot();
	double r = HA.getR();
	vector<double> z1(n), v(n + 1), cmc(N);
	double xi, sigma, d1, d2;
	
	double tmp = 0;
	for (int i = 0; i < N; i++) {
		z1 = rng::rnorm(n);
		v = HA.generateVolPath(T, n, z1);
		v.pop_back(); // remove terminal element
		xi = conditionalXi(rho, v, T/n, z1);
		sigma = conditionalSigma(v, T/n);
		d1 = conditionalD1(S0, K, r, sigma, rho, xi, T);
		d2 = conditionalD2(S0, K, r, sigma, rho, xi, T);
		cmc.at(i) = S0 * xi * MathHelp::stdNormCdf(d1) - 
			K * exp(-r * T) * MathHelp::stdNormCdf(d2);
		tmp += v.size();
	}	
	tmp /= N;
	cout << tmp << endl;
	return cmc;
}
int main() {
	srand(402);	
	double Z_a_2 = 1.96; // two tailed Z value for 1 - alpha confidence level
	// 95% = 1.96, 99% = 2.575, 99.5% = 2.81, 99.9% = 3.29, 99.99% = 3.89
	int log10Nmin = 2;
	int log10Nmax = 6;

	vector<int> N (log10Nmax - log10Nmin + 1); // vector of N to compute
	for (int i = 0; i < N.size(); i++) {
		N.at(i) = pow(10, log10Nmin + i);
	}
			
	// some parameters for the heston model
	double S0 = 100;
	double K = 100;
	char type = 'c';
	double v0 = 0.1;
	double t = 0;
	double T = 0.5;
	int n = 125;
	
	double r = 0;
	double rho_v[] = {-0.7, 0, 0.7};
	string rho_s[] = {"-7","0","7"};
	double kappa = 2;
	double theta = 0.03;
	double sigma = 0.1;
		
	for (int k = 0; k < sizeof(rho_v)/sizeof(rho_v[0]); k++) {
		double rho = rho_v[k];
		string filepath = "../../data/q4_heston_cmc_est";
		filepath = filepath + "_" + rho_s[k] + ".csv";	
		cout << filepath << endl;
	
		// build heston asset w/ given params
		HestonAsset HA(S0, r, v0, sigma, kappa, theta, rho);
		// start conditional sim
		vector<vector<double> > hest_cmc_CI( N.size() ); // confidence intervals over all N
		clock_t t1, t2; float diff;
	
		for (int i = 0; i < N.size(); i++) {
			t1 = clock(); // time code
			vector<double> cmc_price = conditionalMC(HA, K, T, n, N.at(i));
			t2 = clock();
			diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
			cout << N.at(i) << "\t\t" << printf("%.5f", diff) << endl;
			hest_cmc_CI.at(i) = CI( cmc_price, Z_a_2 );
			hest_cmc_CI.at(i).push_back( pow(MathHelp::sd(cmc_price), 2) );
			hest_cmc_CI.at(i).push_back( N.at(i) );
			hest_cmc_CI.at(i).push_back( diff );
		}
		cout << "\n";

		string header = "CI_lo,mean,CI_hi,var,N,time,";
		writeMatrix(hest_cmc_CI, filepath, header);
	}

}





























