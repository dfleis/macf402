#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>
#include "rng.h"
#include "MathHelp.h"
#include "Asset.h"
#include "VanillaOption.h"
#include "AriAsianOption.h"
#include "GeoAsianOption.h"

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
vector<double> crudeSim(VanillaOption VO, double T, int N) {	
	double K = VO.getStrike();
	Asset A = VO.getAsset();
	double r = A.getR();
	vector<double> h(N); // payoff vector
	double z; // random normal variate

	for (int i = 0; i < N; i++) {
		z = rng::rnorm();
		h.at(i) = exp(-r * T) * fmax(K - A.generatePath(T, z), 0);
	}
	return h;
}
double antitheticNormal(double z, double mu) {
	return 2 * mu - z;
}
vector<double> antitheticNormal(vector<double> z, double mu) {
	vector<double> z_AT(z.size());
	for (int i = 0; i < z.size(); i++) {
		z_AT.at(i) = antitheticNormal(z.at(i), mu);
	}
	return z_AT;
}
vector<double> antitheticSim(VanillaOption VO, double T, int N) {
	N /= 2;
	double K = VO.getStrike();
	Asset A = VO.getAsset();
	double r = A.getR();
	vector<double> h(N); // payoff vector
	double z1, z2, S1, S2, h1, h2; // normal variates, asset prices, payoffs
	
	for (int i = 0; i < N; i++) {
		z1 = rng::rnorm();
		z2 = antitheticNormal(z1, 0); // antithetic variate to z1
		S1 = A.generatePath(T, z1);
		S2 = A.generatePath(T, z2);
		h1 = fmax(K - S1, 0);
		h2 = fmax(K - S2, 0);
		h.at(i) = 0.5 * (h1 + h2);
	}
	return h;
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
int main() {
	srand(402);	
	double r = 0.0175; 
	double sigma = 0.4;
	double S0 = 100; 
	double K = 40;
	double t = 0; 
	double T = 1;
	char type = 'p'; 
	int log10Nmin = 2;
	int log10Nmax = 8;
	double Z_a_2 = 1.96; // two tailed Z value for 1 - alpha confidence level
	// 95% = 1.96, 99% = 2.575, 99.5% = 2.81, 99.9% = 3.29, 99.99% = 3.89
	
	Asset A(S0, r, sigma);
	VanillaOption VO(A, K, t, T, type);

	vector<int> N;
	for (int i = log10Nmin; i <= log10Nmax; i++) {
		N.push_back( pow(10, i) );
	}

	vector<vector<double> > crude_CI( N.size() );
 	vector<vector<double> > anti_CI( N.size() );
	clock_t t1, t2; float diff;
	vector<double> h_crude, h_anti;
	
	cout << "crude\t\tantithetic" << endl;
	for (int i = 0; i < N.size(); i++) {
		cout << N.at(i) << "\t\t";
		
		t1 = clock();
		h_crude = crudeSim(VO, T, N.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << "\t\t";
		crude_CI.at(i) = CI(h_crude, Z_a_2);
		crude_CI.at(i).push_back( pow(MathHelp::sd(h_crude), 2) );
		crude_CI.at(i).push_back( N.at(i) );
		crude_CI.at(i).push_back( diff );
		
		t1 = clock();
		h_anti = antitheticSim(VO, T, N.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << endl;
		anti_CI.at(i) = CI(h_anti, Z_a_2);
		anti_CI.at(i).push_back( pow(MathHelp::sd(h_anti), 2) );
		anti_CI.at(i).push_back( N.at(i) );
		anti_CI.at(i).push_back( diff );
	}
 	
	string header = "CI_lo,mean,CI_hi,var,N,time,";
	string filepath_crude = "../../data/q1_crude_est.csv";
 	string filepath_anti = "../../data/q1_anti_est.csv";

	writeMatrix(crude_CI, filepath_crude, header);	
	writeMatrix(anti_CI, filepath_anti, header);
}




















































