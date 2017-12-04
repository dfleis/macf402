/*
	Header file for an AsianOption
*/
#ifndef ASIANOPTION_H
#define ASIANOPTION_H

#include "Option.h"
#include "Asset.h"
#include <vector>

class AsianOption: public Option {
	protected:
		int n;
		std::vector<double> path;
		
	public:
 		AsianOption(); // default constructor
		AsianOption(Asset A, double K, double t, double T, int n, char type);
		double getNMonitoringPoints();
		void setNMonitoringPoints(int n);
		std::vector<double> getPath();
		void setPath(std::vector<double> path);
		
		virtual double payoff() = 0;	
		virtual double price() = 0;	
};
#endif