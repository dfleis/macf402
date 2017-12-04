#include <iostream>
#include <fstream> // file stream
#include <sstream> // string stream
#include <cmath>
#include <ctime> // time code
#include <vector>
#include <numeric> // iota
#include <iomanip> // setfill
using namespace std;

// risk neutral probabilities
const double p = 0.5; // global variable
const double q = 1 - p; // global variable

vector<string> splitstr(string str, char c) { // split string into vector by commas
	vector<string> out;
	stringstream ss(str);
	
	while (ss.good()) {
		string substr;
		getline(ss, substr, c);
		out.push_back(substr);
	}
	return out;
}
double u_fxn(double R, double sigma, double dt) { // up factor function
	return exp(sigma * sqrt(dt) + (R - 0.5 * pow(sigma, 2)) * dt);
}
double d_fxn(double R, double sigma, double dt) { // down factor function
	return exp(-sigma * sqrt(dt) + (R - 0.5 * pow(sigma, 2)) * dt);
}
double price_asset(double S0, int n, int h, double u, double d) { // asset price function
	return S0 * pow(u, h) * pow(d, n - h); // we should probably check if 0 <= h <= n...
}
double payoff(double S, double K, char type) {
	if (type == 'c')  // payoff if call option
		return fmax(S - K, 0); 
	else if (type == 'p') // payoff if put option
		return fmax(K - S, 0);
	else 
		return -1; // do something obviously wrong if input is invalid
}
double price_option_am_r(int n, int h, int N, double S0, double K, 
	double R, double sigma, double dt, char type) { // naive recursive implementation
	// american option price
	
	double u = u_fxn(R, sigma, dt);
	double d = d_fxn(R, sigma, dt);
	double S = price_asset(S0, n, h, u, d);
	double G_n = payoff(S, K, type);
	
	if (n == N) { // terminal node is simply the payoff
		return G_n;
	}
	else { // if we're not at the terminal node, use the recursive algorithm
		double Vu = price_option_am_r(n + 1, h + 1, N, S0, K, R, sigma, dt, type);
		double Vd = price_option_am_r(n + 1, h, N, S0, K, R, sigma, dt, type);
		return fmax( G_n, pow(1 + R, -dt) * ( p * Vu + q * Vd) );
	}
}
double price_option_am_i(int N, double S0, double K, 
	double R, double sigma, double dt, char type) { // iterative implementation
	// american option price
	// assumes you want the time-zero price of the option
	
	double prices[N + 1], S;
	double u = u_fxn(R, sigma, dt);
	double d = d_fxn(R, sigma, dt);
	
	// first loop over all the terminal nodes for all possible H and T combinations
	// to compute V_N = payoff 
	for (int i = 0; i < N + 1; i++) { // traverse upwards through the nodes
		S = price_asset(S0, N, i, u, d); 
		prices[i] = payoff(S, K, type);
	}
	// use the recursive algorithm to price nodes prior to the terminal nodes
	for (int i = N - 1; i >= 0; i--) // traverse backwards through the tree
		for (int j = 0; j < i + 1; j++) { // traverse upwards through the nodes
			
			S = price_asset(S0, i, j, u, d); // price S at i-th step given j heads
			
			prices[j] = fmax( payoff(S, K, type), 
				pow(1 + R, -dt) * (p * prices[j + 1] + q * prices[j]) );
			
		}
	return prices[0];
}

int main() {
	int by = 50;
	int from = 50;
	int to = 3 * pow(10,3);
	int Nvec[to/by]; // number of steps for the binomial model
	for (int i = 0; i < sizeof(Nvec)/sizeof(Nvec[0]); i++) {
		Nvec[i] = from + i * by;
	}
	
	double R = 0.005; // riskless interest rate
	
	string infilepath = "../data/implied_vols_googl_NMAX.csv";
	// INPUT DATA COLUMN HEADERS
	// N, spots, strike, tau, type, ask, vol, time
	
	double dt; char type;
 	vector<double> spots, strikes, taus, asks, vols;
 	vector<string> types;
 	 	
 	// READ DATA
 	/*
 		Our approach is to loop over each line of the CSV and save each column
 		to a corresponding vector. Later we will loop through each vector to compute
		the relevant price
 	*/
	ifstream data( infilepath ); // open csv
	string line; // each row of the csv prior to processing
	vector<string> row; // we will separate the csv rows into their entries
	
	getline(data, line); // read header	
	while (getline(data, line)) { // loop over each line in the CSV
 		row = splitstr(line, ','); // split line on comma
 		 
		spots.push_back( stod(row.at(1)) ); 
		strikes.push_back( stod(row.at(2)) );
		taus.push_back( stod(row.at(3)) );
		types.push_back( row.at(4) );
		asks.push_back( stod(row.at(5)) );
		vols.push_back( stod(row.at(6)) );
	}
  	data.close(); // close input csv

	// COMPUTE AMERICAN OPTION PRICES
	for (int j = 0; j < sizeof(Nvec)/sizeof(Nvec[0]); j++) {
		int N = Nvec[j];
		vector<double> prices, diffs;
		
		for (int i = 0; i < spots.size(); i++) {
			// compute time step value
			dt = taus.at(i) / N;
			
			// extract option type information
			if ( types.at(i) == "C" ) 
				type = 'c';
			else if ( types.at(i) == "P" )
				type = 'p';
			else // do something obviously wrong if input is unexpected
				type = 'z';
	
			// compute price
			clock_t t1 = clock();
			prices.push_back( price_option_am_i(N, spots.at(i), strikes.at(i), 
				R, vols.at(i), dt, type) );	
			clock_t t2 = clock();
			double elapsed = double(t2 - t1) / CLOCKS_PER_SEC;
			cout << "idx: " << i << "\t" << "N = " << N << "\t" << elapsed << endl;
			
			// difference between computed price and observed market price
			diffs.push_back( prices.at(i) - asks.at(i) ); 
		}
	
 		// WRITE DATA
 		stringstream Nss; Nss << setw(4) << setfill('0') << N;
 		string Ns = Nss.str();
		string outfilename = "../data/out/american_prices_googl_N" + Ns + ".csv";
		cout << outfilename << endl;
			
 	 	ofstream outfile;
  		outfile.open( outfilename ); // initialize csv
 	 	string header = "N,spots,strike,tau,type,ask,vol,price,diff";
  		outfile << header << endl; // header
  	  	  	
  		for (int i = 0; i < spots.size(); i++) { // loop over every option in the data set
  			// print to CSV
 	 		outfile << N << "," << spots.at(i) << "," << strikes.at(i) 
  			<< "," << taus.at(i) << "," << types.at(i) << "," << asks.at(i) 
 	 		<< "," << vols.at(i) << "," << prices.at(i) << "," << diffs.at(i) << endl;		
  		}
  	
	  	outfile.close(); // close output csv
	  	cout << "---------------------------------" << endl;
	}
}





































	