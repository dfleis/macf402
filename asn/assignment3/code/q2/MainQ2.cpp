#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "MathHelp.h"
#include "Asset.h"
#include "CallOption.h"
#include "PutOption.h"
#include "OptionData.h"

using namespace std;

vector<double> concatenateVec(vector<double> A, vector<double> B) {
	// concatenate vector A with vector B
	vector<double> AB; 

	AB.reserve( A.size() + B.size() ); // preallocate memory
	AB.insert( AB.end(), A.begin(), A.end() );
	AB.insert( AB.end(), B.begin(), B.end() );
	
	return AB;
}

vector<vector<double> > cartProd(vector<vector<double> > x, vector<vector<double> > y) {
	// cartesian product of two vectors 
	// (well, technically for vectors of vectors for generality)
	int n = x.size();
	int m = y.size();
	vector<vector<double> > out(n + m);
		
	for (int i = 0; i < x.at(0).size(); i++) {
		for (int j = 0; j < y.at(0).size(); j++) {
			for (int k = 0; k < n; k++) {
				out.at(k).push_back( x.at(k).at(i) );
			}
			int idx = 0;
			for (int k = n; k < (n + m); k++) {
				out.at(k).push_back( y.at(idx).at(j) );
			}
		}
	}
	return out;
}
vector<vector<double> > cartProdN(vector<vector<double> > data) {
	// generalizes to the cartesian product over N vectors
	int d = data.size(); // d dimensions/vectors to perform the cartesian product
	
	if (d == 1) { // you're killing me here, come on
		return data;
	}
	else { // if the number of dimensions/vectors >= 2
		vector<vector<double> > x1(1);
		vector<vector<double> > x2(1);
		vector<vector<double> > xn(1);

		x1.at(0) = data.at(0);
		x2.at(0) = data.at(1);
		
		vector<vector<double> > data_cart_prod = cartProd(x1, x2);
		
		for (int i = 2; i < d; i++) {
			xn.at(0) = data.at(i);
			data_cart_prod = cartProd(data_cart_prod, xn);
		}
		return data_cart_prod;
	}
}

vector<vector<double> > fitVols(double spot, vector<vector<double> > x_obs, 
	vector<double> y_obs, vector<double> h_vec, vector<double> min_vec,
	vector<double> max_vec, int N) { 
	// h_vec, min_vec, max_vec, and number of x_obs dimensions must agree
	
	int d = h_vec.size();
	vector<double> dx(d); // granularity of each x dimension to fit
	vector<vector<double> > x_unobs(d); // matrix of points to smooth over
	// ideally i would have initialized the dimension of each 'column' in the matrix
	// using N_vec, but i don't know how to do so (thus forced to use push_back)
			
	for (int i = 0; i < d; i++) { // compute dx for each x var
		dx.at(i) = (max_vec.at(i) - min_vec.at(i))/(N - 1);
		x_unobs.at(i).push_back( min_vec.at(i) ); // place the min as the first step
	}
		
	for (int i = 0; i < d; i++) { // generate matrix of steps, loop over each x var
		for (int j = 1; j < N; j++) { // increment over corresp. dx
			x_unobs.at(i).push_back( x_unobs.at(i).at(j - 1) + dx.at(i) ); 
		}
	}
		
	// create the cartesian product of our x coordinates to be smoothed over
	vector<vector<double> > x_unobs_cart_prod = cartProdN(x_unobs);
	
	vector<vector<double> > data_out;
	vector<double> row(d + 1); // create vector for each 'row' of the data 
	vector<double> x_unobs_row(d); // single row of our cartesian product
	
	for (int i = 0; i < x_unobs_cart_prod.at(0).size(); i++) {
		for (int j = 0; j < x_unobs_cart_prod.size(); j++) {
			row.at(j) = x_unobs_cart_prod.at(j).at(i);
			x_unobs_row.at(j) = x_unobs_cart_prod.at(j).at(i);
		}
		row.at(d) = MathHelp::nadWatEstimator(h_vec, x_unobs_row, x_obs, y_obs);
		
		data_out.push_back( row ); // append row to data output
	}
	return data_out;
}

vector<vector<double> > fitVols(double spot, vector<double> x_obs, vector<double> y_obs,
	double h, double min, double max, int N) {
	// 1 dimensional wrapper
	// we are storing our 1-D parameters as vectors since we have defined
	// our smoothing fxns to be sufficiently general
	
	vector<double> h_vec(1); h_vec.at(0) = h; // smoothing parameter
	vector<double> min_vec(1), max_vec(1);
	min_vec.at(0) = min; max_vec.at(0) = max; // min/max x-variable to fit
	

	// independent variables: place observed x variables in matrix
	vector<vector<double> > x_obs_vec(1);
	x_obs_vec.at(0) = x_obs;
	
	vector<vector<double> > fitted_vols = fitVols(spot, x_obs_vec, 
		y_obs, h_vec, min_vec, max_vec, N);
	
	return fitted_vols;
}

int main() {
	// ----------------------
	// set up some parameters
	// ----------------------
	
	double spot = 502.12; // market spot
	string start_date = "2012-02-17"; // t0
	double r = 0.0175; // con't comp. interest rate
	
	Asset AAPL(spot, r); // create our asset
		
	double bisect_a = 0; // bisection left endpoint
	double bisect_b = 2; // bisection right endpoint
	
	// ---------
	// load data
	// ---------
	
	// filenames to load data
	string filename_mar2012 = "../../data/cleaned/cleaned_mar2012.csv";
	string filename_apr2012 = "../../data/cleaned/cleaned_apr2012.csv";
	string filename_may2012 = "../../data/cleaned/cleaned_may2012.csv";
		
	// load/create datasets
	OptionData mar2012, apr2012, may2012;
	mar2012.readCsv(filename_mar2012); 
	apr2012.readCsv(filename_apr2012); 
	may2012.readCsv(filename_may2012);
	
	// -----------
	// filter data
	// -----------
	
	// keep only data with volume >= 100, last >= 0.1
	double inf = numeric_limits<double>::infinity();
	mar2012.filterLasts(0.1, inf); mar2012.filterVolumes(100, inf);
	apr2012.filterLasts(0.1, inf); apr2012.filterVolumes(100, inf);
	may2012.filterLasts(0.1, inf); may2012.filterVolumes(100, inf);
	
	// generate maturities & implied vols
	mar2012.computeTaus(start_date);
	apr2012.computeTaus(start_date);
	may2012.computeTaus(start_date);
	mar2012.computeImpliedVols(AAPL, bisect_a, bisect_b);
	apr2012.computeImpliedVols(AAPL, bisect_a, bisect_b);
	may2012.computeImpliedVols(AAPL, bisect_a, bisect_b);
	
	// ---------------------------------------
	// write unsmoothed data to csv
	// ---------------------------------------
	
	// prepare csv files
	ofstream mar_outfile, apr_outfile, may_outfile;
	mar_outfile.open("../../data/output/q2_imp_vols_mar.csv"); // initialize csv
	apr_outfile.open("../../data/output/q2_imp_vols_apr.csv"); // initialize csv
	may_outfile.open("../../data/output/q2_imp_vols_may.csv"); // initialize csv

	// write csv columns & write data
	mar_outfile << "date,maturity,last,strike,moneyness,type,imp_vol" << endl;
	apr_outfile << "date,maturity,last,strike,moneyness,type,imp_vol" << endl;
	may_outfile << "date,maturity,last,strike,moneyness,type,imp_vol" << endl;

	// loop over march 2012 rows
	for (int i = 0; i < mar2012.getLasts().size(); i++) {
		string date = mar2012.getDates().at(i);
		double tau = mar2012.getTaus().at(i);
		double last = mar2012.getLasts().at(i);
		double vol = mar2012.getVolumes().at(i);
		double strike = mar2012.getStrikes().at(i);
		string type = mar2012.getTypes().at(i);
		double imp_vol = mar2012.getImpliedVols().at(i);
		
		mar_outfile << date << "," << tau << "," << last << "," << strike << "," << 
			spot/strike << "," << type << "," << imp_vol << endl;	
	}

	// loop over april 2012 rows
	for (int i = 0; i < apr2012.getLasts().size(); i++) {
		string date = apr2012.getDates().at(i);
		double tau = apr2012.getTaus().at(i);
		double last = apr2012.getLasts().at(i);
		double vol = apr2012.getVolumes().at(i);
		double strike = apr2012.getStrikes().at(i);
		string type = apr2012.getTypes().at(i);
		double imp_vol = apr2012.getImpliedVols().at(i);
		
		apr_outfile << date << "," << tau << "," << last << "," << strike << "," << 
			spot/strike << "," << type << "," << imp_vol << endl;	
	}
	
	// loop over may 2012 rows
	for (int i = 0; i < may2012.getLasts().size(); i++) {
		string date = may2012.getDates().at(i);
		double tau = may2012.getTaus().at(i);
		double last = may2012.getLasts().at(i);
		double vol = may2012.getVolumes().at(i);
		double strike = may2012.getStrikes().at(i);
		string type = may2012.getTypes().at(i);
		double imp_vol = may2012.getImpliedVols().at(i);
		
		may_outfile << date << "," << tau << "," << last << "," << strike << "," << 
			spot/strike << "," << type << "," << imp_vol << endl;	
	}
	
	mar_outfile.close(); apr_outfile.close(); may_outfile.close();	
	
	// -------------
	// 1-D Smoothing
	// -------------
	double h_1d = 0.075; // smoothing parameter
	double min_1d = 0.5; // min moneyness to fit
	double max_1d = 1.5; // max moneyness to fit
	int N_1d = 200; // granularity of moneyness steps to fit
	
	// --------------
	// Fit March Data
	// --------------
	// prepare march moneyness
	vector<double> moneyness_mar(mar2012.getStrikes().size());
	for (int i = 0; i < moneyness_mar.size(); i++) {
		moneyness_mar.at(i) = spot / mar2012.getStrikes().at(i);
	}
	
	// fit vols
	vector<vector<double> > fitted_vols_1d_mar = fitVols(spot, moneyness_mar, 
		mar2012.getImpliedVols(), h_1d, min_1d, max_1d, N_1d);

	// --------------
	// Fit April Data
	// --------------	// prepare april moneyness
	vector<double> moneyness_apr(apr2012.getStrikes().size());
	for (int i = 0; i < moneyness_apr.size(); i++) {
		moneyness_apr.at(i) = spot / apr2012.getStrikes().at(i);
	}
	
	// fit vols
	vector<vector<double> > fitted_vols_1d_apr = fitVols(spot, moneyness_apr, 
		apr2012.getImpliedVols(), h_1d, min_1d, max_1d, N_1d);

	// --------------
	// Fit May Data
	// --------------	// prepare may moneyness
	vector<double> moneyness_may(may2012.getStrikes().size());
	for (int i = 0; i < moneyness_may.size(); i++) {
		moneyness_may.at(i) = spot / may2012.getStrikes().at(i);
	}
	
	// fit vols
	vector<vector<double> > fitted_vols_1d_may = fitVols(spot, moneyness_may, 
		may2012.getImpliedVols(), h_1d, min_1d, max_1d, N_1d);

	// -------------
	// Write to CSV
	// -------------
	// prepare csv
	ofstream fitted_1d_mar_outfile, fitted_1d_apr_outfile, fitted_1d_may_outfile;
	fitted_1d_mar_outfile.open("../../data/output/q2_fitted_1d_mar_data.csv");
	fitted_1d_apr_outfile.open("../../data/output/q2_fitted_1d_apr_data.csv");
	fitted_1d_may_outfile.open("../../data/output/q2_fitted_1d_may_data.csv");

	// write headers
	fitted_1d_mar_outfile << "moneyness,fitted_vol" << endl;
	fitted_1d_apr_outfile << "moneyness,fitted_vol" << endl;
	fitted_1d_may_outfile << "moneyness,fitted_vol" << endl;

	// write rows
 	for (int i = 0; i < fitted_vols_1d_mar.size(); i++) { // loop over march rows
 		fitted_1d_mar_outfile << fitted_vols_1d_mar.at(i).at(0) << "," <<
  			fitted_vols_1d_mar.at(i).at(1) << endl;
 	}
	for (int i = 0; i < fitted_vols_1d_apr.size(); i++) { // loop over april rows
 		fitted_1d_apr_outfile << fitted_vols_1d_apr.at(i).at(0) << "," <<
  			fitted_vols_1d_apr.at(i).at(1) << endl;
 	}
 	for (int i = 0; i < fitted_vols_1d_may.size(); i++) { // loop over may rows
 		fitted_1d_may_outfile << fitted_vols_1d_may.at(i).at(0) << "," <<
  			fitted_vols_1d_may.at(i).at(1) << endl;
 	}
 	
 	// close files
 	fitted_1d_mar_outfile.close();
	fitted_1d_apr_outfile.close();
	fitted_1d_may_outfile.close();
	
	// -------------
	// 2-D Smoothing
	// -------------	
	vector<double> h_2d(2); // smoothing parameters
	h_2d.at(0) = 0.05; 
	h_2d.at(1) = 0.05;
	vector<double> min_2d(2), max_2d(2);
	min_2d.at(0) = spot/770; min_2d.at(1) = 29.0/365; // min moneyness/time to fit
	max_2d.at(0) = spot/330; max_2d.at(1) = 90.0/365; // max moneyness/time to fit
	int N_2d = 50;
	
	// concatenate our vectors of strikes for use in the smoother
	vector<double> mar_apr_strikes = concatenateVec(mar2012.getStrikes(), 
		apr2012.getStrikes());
	vector<double> strikes = concatenateVec(mar_apr_strikes, may2012.getStrikes());
	
	// compute moneyness 
	vector<double> moneyness(strikes.size());
	for (int i = 0; i < strikes.size(); i++) {
		moneyness.at(i) = spot/strikes.at(i);
	}
	
	// concatenate our vectors of maturities for use in the smoother
	vector<double> mar_apr_taus = concatenateVec(mar2012.getTaus(), 
		apr2012.getTaus());
	vector<double> taus = concatenateVec(mar_apr_taus, may2012.getTaus());
	
	// independent variable: place all x variables in a matrix
	vector<vector<double> > x_obs(2);
	x_obs.at(0) = moneyness;
	x_obs.at(1) = taus;
	
	// dependent variable: create vector of implied vols
	vector<double> mar_apr_vols = concatenateVec(mar2012.getImpliedVols(), 
		apr2012.getImpliedVols());
	vector<double> vols = concatenateVec(mar_apr_vols, may2012.getImpliedVols());	

	// --------
	// Fit Data
 	// --------
	vector<vector<double> > fitted_vols_2d = fitVols(spot, x_obs, vols, 
		h_2d, min_2d, max_2d, N_2d);

	// -------------
	// Write to CSV
	// -------------
	// prepare to write 2d fitted data to csv
	ofstream fitted_2d_outfile;
	fitted_2d_outfile.open("../../data/output/q2_fitted_2d_data.csv");
	fitted_2d_outfile << "moneyness,maturity,fitted_vol" << endl;
	
	// write rows to csv
 	for (int i = 0; i < fitted_vols_2d.size(); i++) {
 		fitted_2d_outfile << fitted_vols_2d.at(i).at(0) << "," <<
  			fitted_vols_2d.at(i).at(1) << "," <<
  			fitted_vols_2d.at(i).at(2) << endl;
 	}
 	
 	fitted_2d_outfile.close();
}






































