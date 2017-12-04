/*
	Header file for a OptionData structure
	Very rudimentary data type to store our option data
*/
#ifndef OPTIONDATA_H
#define OPTIONDATA_H

#include "Option.h"
#include <ctime> // calculate time between dates
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "Asset.h"

using namespace std;

class OptionData {
	private:
		vector<string> dates;
		vector<double> lasts;
		vector<double> volumes;
		vector<double> strikes;
		vector<string> types;
		
		vector<double> taus;
		vector<double> implied_vols;
		
		vector<int> stringToDate(string string_date);
		tm makeTm(vector<int> date);
		double stringDateDiff(string date_a, string date_b);

	public:
 		OptionData(); // default constructor
		OptionData(vector<string> dates_in, vector<double> lasts_in, 
			vector<double> volumes_in, vector<double> strikes_in, vector<string> types_in);
		vector<string> getDates();
		vector<double> getLasts(); 
		vector<double> getVolumes();
		vector<double> getStrikes();
		vector<string> getTypes();
		vector<double> getTaus();
		vector<double> getImpliedVols();
					
		void readCsv(string filename);
		void filterVolumes(double min, double max);
		void filterLasts(double min, double max);
		
		void computeTaus(string start_date);
		void computeImpliedVols(Asset A, double bisect_a, double bisect_b);
};
#endif