#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime> // calculate time between dates
#include <time.h> // time code
#include "Asset.h"
#include "CallOption.h"
#include "PutOption.h"
#include "MathHelp.h"

int main() {	
	double bisect_a = 0; // left endpoint for starting bisection
	double bisect_b = 1; // right endpoint for starting bisection 
	// double NR_x0 = 0.1; // initial guess for newton-raphson
	
	double S0 = 4.75; // market spot price
	double R = 0.0492; // con't comp. interest rate
	
	// specify the option chain strikes & asks
	double strikes[8] = {4,4.25,4.5,4.75,5,5.25,5.5,5.75};
	double call_asks[8] = {nan("1"),nan("1"),0.33,0.16,0.06,0.02,0.01,0.01};
	double put_asks[8] = {0.02,0.04,0.09,0.20,0.38,0.59,nan("1"),nan("1")};
		
	// compute tau from date endpoints
	struct std::tm date_a = {0,0,0,29,6,102}; // July 29th 2002
	struct std::tm date_b = {0,0,0,26,8,102}; // Sept 26th 2002
	std::time_t start_date = std::mktime(&date_a);
	std::time_t end_date= std::mktime(&date_b);
	double tau = std::difftime(end_date, start_date)/(60 * 60 * 24 * 365);
	
	Asset A(S0, R); // create our underlying asset
	
	// prepare to write data to csv
	std::ofstream outfile_vols;
	outfile_vols.open("../../data/output/q1_implied_vols.csv"); // initialize csv

	// write csv columns & write data
	outfile_vols << "strike,moneyness,call_imp_vol,put_imp_vol" << std::endl;
	
	// loop over strike & ask values & compute implied volatility estimates
	for (int i = 0; i < sizeof(strikes)/sizeof(strikes[0]); i++) {
		
		CallOption C(A, strikes[i], tau);
		PutOption P(A, strikes[i], tau);
			
		CallOption* CO = &C;
		PutOption* PO = &P;
						
		// Newton-Raphson
		// outfile_vols << strikes[i] << "," << strikes[i]/S0  << "," <<
		//	CO -> implied_vol_NR(call_asks[i], NR_x0) << "," <<
		// 	PO -> implied_vol_NR(put_asks[i], NR_x0) << std::endl;
		
		// Bisection
		outfile_vols << strikes[i] << "," << strikes[i]/S0  << "," <<
			CO -> impliedVolBisect(call_asks[i], bisect_a, bisect_b) << "," <<
			PO -> impliedVolBisect(put_asks[i], bisect_a, bisect_b) << std::endl;
	}
		
	outfile_vols.close();
}


















	