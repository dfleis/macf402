#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>
#include "rng.h"
#include "MathHelp.h"
#include "Asset.h"
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
vector<double> crudeSim(AsianOption* AO, int N) {
	Asset A = AO->getAsset();
	int n = AO->getNMonitoringPoints();
	double delta_t = AO->getTEnd() - AO->getTStart();
	
	vector<double> h(N), z(n);
	
	for (int i = 0; i < N; i++) {
		z = rng::rnorm(n);
		AO->setPath( A.generatePath(delta_t, n, z) );
		h.at(i) = AO->payoff();
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
vector<double> antitheticSim(AsianOption* AO, int N) {
	N /= 2;
	Asset A = AO->getAsset();
	int n = AO->getNMonitoringPoints();
	double delta_t = AO->getTEnd() - AO->getTStart();
	
	vector<double> h(N), z1(n), z2(n);
	double h1, h2;
	
	for (int i = 0; i < N; i++) {
		z1 = rng::rnorm(n);
		z2 = antitheticNormal(z1, 0);		
		AO->setPath( A.generatePath(delta_t, n, z1) );
		h1 = AO->payoff();
		AO->setPath( A.generatePath(delta_t, n, z2) );
		h2 = AO->payoff();
		h.at(i) = 0.5 * (h1 + h2);
	}
	return h;
}
double controlC(std::vector<double> x, std::vector<double> y) {
	int m = x.size(); // we should check if x, y dimensions agree
	double theta_hat_x = MathHelp::mean(x);
	double theta_hat_y = MathHelp::mean(y);
	double s_hat_y = MathHelp::sd(y);

	double z = 0;
	for (int i = 0; i < m; i++) {
		z += x.at(i) * y.at(i);
	}
	return (z - m * theta_hat_x * theta_hat_y) / ((m - 1) * pow(s_hat_y, 2));
}
vector<double> controlSim(AriAsianOption AAO, GeoAsianOption GAO, int N, int m) {
	int M = N - m;
	Asset A = AAO.getAsset(); // check if assets are the same
	int n = AAO.getNMonitoringPoints();
	double delta_t = AAO.getTEnd() - AAO.getTStart();
	
	vector<double> z(n), S(n + 1), p1(m), p2(m), h(M);
	for (int i = 0; i < m; i++) { // M = nb of pilot runs to compute c
		z = rng::rnorm(n);
		S = A.generatePath(delta_t, n, z);
		AAO.setPath(S);
		GAO.setPath(S);
		p1.at(i) = AAO.payoff();
		p2.at(i) = GAO.payoff();
	}
	double c = controlC(p1, p2);

	for (int i = 0; i < M; i++) {
		z = rng::rnorm(n);
		S = A.generatePath(delta_t, n, z);
		AAO.setPath(S);
		GAO.setPath(S);
		h.at(i) = AAO.payoff() - c * ( GAO.payoff() - GAO.price() );
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
	double sigma = 0.25;
	double S0 = 500; 
	double K = 500;
	double t = 0; 
	double T = 1;
	int n = 52; 
	char type = 'c'; 
	int log10Nmin = 2;
	int log10Nmax = 8;
	double Z_a_2 = 1.96; // two tailed Z value for 1 - alpha confidence level
	// 95% = 1.96, 99% = 2.575, 99.5% = 2.81, 99.9% = 3.29, 99.99% = 3.89
	
	Asset A(S0, r, sigma);
	AriAsianOption AAO(A, K, t, T, n , type);	
	GeoAsianOption GAO(A, K, t, T, n, type);
	double geo_price = GAO.price();
		
	vector<int> N, m;
	for (int i = log10Nmin; i <= log10Nmax; i++) {
		N.push_back( pow(10, i) );
		m.push_back( pow(10, i - 1) );
	}
		
	vector<vector<double> > crude_CI( N.size() );
	vector<vector<double> > anti_CI( N.size() );
	vector<vector<double> > cont_CI( N.size() );
	vector<vector<double> > geo_crude_CI( N.size() );
	vector<vector<double> > geo_anti_CI( N.size() );
	
	clock_t t1, t2; float diff;
	vector<double> h_crude, h_anti, h_cont, geo_h_crude, geo_h_anti;
	cout << "crude\tgeo_crude\tanti\tgeo_anti\tcontrol" << endl;
	for (int i = 0; i < N.size(); i++) {
		t1 = clock();
		h_crude = crudeSim(&AAO, N.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << "\t";
		crude_CI.at(i) = CI(h_crude, Z_a_2);
		crude_CI.at(i).push_back( pow(MathHelp::sd(h_crude), 2) );
		crude_CI.at(i).push_back( geo_price );
		crude_CI.at(i).push_back( N.at(i) );
		crude_CI.at(i).push_back( diff );
		
		t1 = clock();
		geo_h_crude = crudeSim(&GAO, N.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << "\t";
		geo_crude_CI.at(i) = CI(geo_h_crude, Z_a_2);
		geo_crude_CI.at(i).push_back( pow(MathHelp::sd(geo_h_crude), 2) );
		geo_crude_CI.at(i).push_back( geo_price );
		geo_crude_CI.at(i).push_back( N.at(i) );
		geo_crude_CI.at(i).push_back( diff );
		
		t1 = clock();
		h_anti = antitheticSim(&AAO, N.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << "\t";
		anti_CI.at(i) = CI(h_anti, Z_a_2);
		anti_CI.at(i).push_back( pow(MathHelp::sd(h_anti), 2) );
		anti_CI.at(i).push_back( geo_price );
		anti_CI.at(i).push_back( N.at(i) );
		anti_CI.at(i).push_back( diff );
		
		t1 = clock();
		geo_h_anti = antitheticSim(&GAO, N.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << "\t";
		geo_anti_CI.at(i) = CI(geo_h_anti, Z_a_2);
		geo_anti_CI.at(i).push_back( pow(MathHelp::sd(geo_h_anti), 2) );
		geo_anti_CI.at(i).push_back( geo_price );
		geo_anti_CI.at(i).push_back( N.at(i) );
		geo_anti_CI.at(i).push_back( diff );

		t1 = clock();
		h_cont = controlSim(AAO, GAO, N.at(i), m.at(i));
		t2 = clock();
		diff = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
		cout << printf("%.5f", diff) << endl;
		cont_CI.at(i) = CI(h_cont, Z_a_2);
		cont_CI.at(i).push_back( pow(MathHelp::sd(h_cont), 2) );
		cont_CI.at(i).push_back( geo_price );
		cont_CI.at(i).push_back( N.at(i) );
		cont_CI.at(i).push_back( diff );
	}
	
	string header = "CI_lo,mean,CI_hi,var,geo_price,N,time,";
	string filepath_crude = "../../data/q3_crude_est.csv";
	string filepath_anti = "../../data/q3_anti_est.csv";
	string filepath_cont = "../../data/q3_cont_est.csv";
	string filepath_geo_crude = "../../data/q3_geo_crude_est.csv";
	string filepath_geo_anti = "../../data/q3_geo_anti_est.csv";

	writeMatrix(crude_CI, filepath_crude, header);	
	writeMatrix(anti_CI, filepath_anti, header);
	writeMatrix(cont_CI, filepath_cont, header);
	writeMatrix(geo_crude_CI, filepath_geo_crude, header);	
	writeMatrix(geo_anti_CI, filepath_geo_anti, header);
}




















































