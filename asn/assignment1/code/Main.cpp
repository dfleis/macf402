#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime> // calculate time between dates
#include <time.h> // time code
#include "Asset.h"
#include "BinomialModel.h"
#include "CallOption.h"
#include "PutOption.h"
#include "MathHelp.h"

int main() {
// 	double ask = 0.33;
// 	double S0 = 4.75;
// 	double strike = 4.5;
// 	double R = 0.0492;
// 	double tau = 0.161;
// 	int N = 100;
// 	
// 	Asset A(S0, R);
// 	CallOption C(A, strike, tau);
// 	CallOption* CO = &C;
// 	
// 	std::cout << BinomialModel::impliedVol(ask, CO, N) << std::endl;
// 	//std::cout << BinomialModel::optionValuation(CO, N) << std::endl;
// 	//std::cout << BinomialModel::assetPrice(A, tau, 5, N) << std::endl;

	double S0 = 4.75; // market spot price
	double R = 0.0492; // con't comp. interest rate
	// steps for the binomial tree calculations
	int N[] = {pow(10,1),pow(10,2),pow(10,3),pow(10,4),pow(10,5),pow(10,6)};
	
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
	std::ofstream outfile_vols, outfile_times;	
	outfile_vols.open("../data/implied_vols.csv"); // initialize csv
	outfile_times.open("../data/times.csv"); // initialize csv

	// write csv columns & write data
	outfile_vols << "N,strike,moneyness,call,put" << std::endl;
	
	// display computation time for each N while looping
	std::cout << "N\t\tseconds" << std::endl;
	outfile_times << "N,seconds" << std::endl;
	
	for (int j = 0; j < (sizeof(N)/sizeof(N[0])); j++) {  // loop over values of N

		clock_t t1, t2;
		t1 = clock();

		// loop over strike & ask values & compute implied volatility estimates
		for (int i = 0; i < sizeof(strikes)/sizeof(strikes[0]); i++) {
		
			CallOption C(A, strikes[i], tau);
			PutOption P(A, strikes[i], tau);
			
			CallOption* CO = &C;
			PutOption* PO = &P;
						
			outfile_vols << N[j] << "," << strikes[i] << "," << strikes[i]/S0  << ","
			<< BinomialModel::impliedVol(call_asks[i],CO,N[j]) << ","
			<< BinomialModel::impliedVol(put_asks[i],PO,N[j]) << std::endl;
		}
		
		t2 = clock();
		float diff ((float)t2 - (float)t1);
		std::cout << "10^" << log10(N[j]) << "\t\t" << diff/CLOCKS_PER_SEC << std::endl;
		outfile_times << N[j] << "," << diff/CLOCKS_PER_SEC << std::endl;
	}
	
	outfile_vols.close();
	outfile_times.close();	
}



	