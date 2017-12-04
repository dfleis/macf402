#include "BinomialModel.h"

double BinomialModel::p = 0.5;
double BinomialModel::u = 1;
double BinomialModel::d = 1;

double BinomialModel::assetPrice(Asset A, double tau, int n, int N) {
	// compute future asset price at time tau given
	// n up steps in a N steps in a binomial tree in tau time
	if (N <= n) return A.getSpot();
	
	double sigma = A.getSigma();
	double R = A.getR();
	double dt = tau/N; // time step interval
	
	// compute up & down factors
	u = exp(sigma * sqrt(dt) + (R - 0.5 * pow(sigma, 2)) * dt);
	d = exp(-sigma * sqrt(dt) + (R - 0.5 * pow(sigma, 2)) * dt);
	
	// Sn = S0 * u^n * d^(N-n)
	double S_N = A.getSpot() * pow(u, n) * pow(d, N - n);
	return S_N;
} 
double BinomialModel::optionValuation(Option* O, int N) {
	// compute option value at time 0 given N steps in a binomial tree
	
	// we rely heavily on pointers since Option O will either contain s
	// CallOption or a PutOption and we wish to use their function implementations
	// so that the code may remain general as possible 
	// (also, I'm also not very good at C++)
	
	Asset A = O->getAsset();
	double R = A.getR();
	double tau = O->getTau();
	double strike = O->getStrike();

	double S_N, V_N;
	double V_0 = 0;
	
	if (N <= 1000) { // if N < 1000 we can use an exact algorithm for the binomial coef.
		for (int k = 0; k <= N; k++) {
			// compute asset value at time N given k up steps
			S_N = assetPrice(A, tau, k, N);
			V_N = O->payoff(S_N);
			
			V_0 += MathHelp::bin_coef(N, k) * pow(p, k) * pow(1 - p, N - k) * V_N;			
		}
	}
	else { // if N too large it then we get into data type representation issues
		   // so we've got to do some trickery
		for (int k = 0; k <= N; k++) {
			// compute asset value at time N given k up steps
			S_N = assetPrice(A, tau, k, N);
			V_N = O->payoff(S_N); // option payoff

			// compute partial sum of the value V_0 of our option
			// exploit the fact that x*y*z = exp(log(x) + log(y)) * z
			// to avoid representation issues
			V_0 += exp(MathHelp::bin_coef(N, k) + k * log(p) + (N - k) * log (1 - p)) * V_N;
		}
	}
	return exp(-R * tau) * V_0; // multiply by discount factor
}
double BinomialModel::impliedVol(double V_obs, Option* O, int N, double eps) {
	// bisection algorithm to back out implied volatility estimate
	// given observed market prices V_obs, Option O, N steps in a binomial model
	// and some tolerance parameter eps
	
	if (V_obs != V_obs) { // exploit the fact that nan is not equal to itself
		// this exists just for simplicity when handling options without asks
		return nan("1");
	}
	
	Asset A = O->getAsset();
	double spot = A.getSpot();
	double R = A.getR();
	
	double a = 0; // vol can never be < 0
	double b = 1; // what is a reasonable upper limit?
	double sigma_guess = (a + b)/2; // bisect!
	
	O->setAsset(Asset(spot, R, sigma_guess));
	double V_guess = optionValuation(O, N);
	
	while (std::abs(V_guess - V_obs) >= eps) {	
		if (V_guess > V_obs) { // guessed price too high, guessed vol too high
			b = sigma_guess; // shrink right endpoint
			sigma_guess = (a + b)/2; // bisect!
		}
		else { // guessed price too low, guessed vol too low
			a = sigma_guess; // shrink left endpoint
			sigma_guess = (a + b)/2; // bisect!
		}
		
		O->setAsset(Asset(spot, R, sigma_guess)); // change value of sigma in underlying asset
		V_guess = optionValuation(O, N); // recompute our guess
	}	
	return sigma_guess;
}

















