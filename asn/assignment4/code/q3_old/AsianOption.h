/*
	Header file for an AsianOption
*/
#ifndef ASIANOPTION_H
#define ASIANOPTION_H

#include "Asset.h"
#include <vector>

class AsianOption {
	protected:
		Asset A;
		double K;
		double t;
		double T;
		int n;
		char type;
		std::vector<double> path;
		
	public:
 		AsianOption(); // default constructor
		AsianOption(Asset A, double K, double t, double T, int n, char type);
		Asset getAsset();
		void setAsset(Asset A);
		double getTStart();
		void setTStart(double t);
		double getTEnd();
		void setTEnd(double T);
		double getNMonitoringPoints();
		void setNMonitoringPoints(int n);
		std::vector<double> getPath();
		void setPath(std::vector<double> path);
		
		virtual double payoff() = 0;	
		virtual double price() = 0;	
};
#endif