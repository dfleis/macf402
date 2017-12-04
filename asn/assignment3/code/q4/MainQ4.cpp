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

vector<vector<double> > readFittedVols(string filename) {
	string fileline;
	vector<vector<double> > data_out(3);

	ifstream data;
	data.open(filename);
	
	if (!data.is_open()) { // abort program if file doesn't exist
		exit(EXIT_FAILURE); // EXIT_FAILRuE comes from <cstdlib>
	}
		
	data >> fileline; // read first line (first line is header data)
	data >> fileline; // read next line
	string last_line;
	vector<string> tmp_line; // temporary storage for each line	
	
	// ASSUMES NO CONSECUTIVE LINES OUGHT TO BE INDENTICAL
	// this solves an issue where the final line was being read twice
	while (!data.eof() || (fileline != last_line)) { // loop over each line of the csv		
		last_line = fileline;

		boost::split(tmp_line, fileline, boost::is_any_of(",")); // split line by comma
		
		data_out.at(0).push_back( stod(tmp_line.at(0)) );
		data_out.at(1).push_back( stod(tmp_line.at(1)) );
		data_out.at(2).push_back( stod(tmp_line.at(2)) );
		
		data >> fileline; // read next line
	}
	return data_out;
}

int main() {
	// ----------------------
	// set up some parameters
	// ----------------------
	
	double spot = 502.12; // market spot
	double r = 0.0175; // con't comp. interest rate
	double strike_wanted = 572.50;
	double maturity_wanted = 60.0/365;
	double moneyness_wanted = spot/strike_wanted;
		
	// ---------
	// load data
	// ---------

	// filenames to load data
	string fitted_vols_infile = "../../data/output/q2_fitted_2d_data.csv";
	
	// load data
	vector<vector<double> > fitted_vols = readFittedVols(fitted_vols_infile);
	
	// --------------------------------
	// find closest fit to desired data
	// --------------------------------
	double closest_moneyness, closest_maturity, sigma_hat, 
		tmp_moneyness_err, tmp_maturity_err;
	
	double moneyness_err = numeric_limits<double>::infinity();
	double maturity_err = numeric_limits<double>::infinity();
	
	for (int i = 0; i < fitted_vols.at(0).size(); i++) { // loop over all values to find
		// the one with the lowest error in both dimensions
		tmp_moneyness_err = fabs(fitted_vols.at(0).at(i) - moneyness_wanted);
		tmp_maturity_err = fabs(fitted_vols.at(1).at(i) - maturity_wanted);
		
		if ( (moneyness_err >=  tmp_moneyness_err) & (maturity_err >= tmp_maturity_err) ) {
			moneyness_err = tmp_moneyness_err;
			maturity_err = tmp_maturity_err;
				
			closest_moneyness = fitted_vols.at(0).at(i);
			closest_maturity = fitted_vols.at(1).at(i);
			
			sigma_hat = fitted_vols.at(2).at(i);
		}
	}
	
	// ------------
	// price option
	// ------------
	Asset AAPL(spot, r, sigma_hat); // create our asset with the sigma_hat
	CallOption C(AAPL, strike_wanted, maturity_wanted); // create our option
	CallOption* CO = &C; // create a pointer to the option
	
	double bs_price = CO -> blackScholesPrice();
	double bs_hedge_ratio = CO -> blackScholesDelta();
	
	cout << closest_moneyness << "\t\t" << closest_maturity << "\t\t" << sigma_hat << endl;
	
	cout << bs_price << endl; // price
	cout << bs_hedge_ratio << endl; // hedge ratio
}






























