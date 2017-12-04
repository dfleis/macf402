#include "OptionData.h"
#include "CallOption.h"
#include "PutOption.h"
#include <iostream>
#include <fstream>

using namespace std;

OptionData::OptionData() {}
OptionData::OptionData(vector<string> dates_in, vector<double> lasts_in,
	vector<double> volumes_in, vector<double> strikes_in, vector<string> types_in) {

	dates = dates_in;
	lasts = lasts_in;
	volumes = volumes_in;
	strikes = strikes_in;
	types = types_in;
}
vector<std::string> OptionData::getDates() {
	return dates;
}
vector<double> OptionData::getLasts() {
	return lasts;
}
vector<double> OptionData::getVolumes() {
	return volumes;		
}
vector<double> OptionData::getStrikes() {
	return strikes;
}
vector<string> OptionData::getTypes() {
	return types;
}
vector<double> OptionData::getTaus() {
	return taus;
}
vector<double> OptionData::getImpliedVols() {
	return implied_vols;
}
			
void OptionData::readCsv(string filename) {
	string fileline;
	
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
			
		dates.push_back(tmp_line[0]); // place date in its vector
		
		// some cells are NA for options that do not exist
		// we must handle these cases
		if (tmp_line[1] != "NA") { // place last in its vector
			lasts.push_back(stod(tmp_line[1]));
		} else {
			lasts.push_back(nan("1"));
		}
		
		if (tmp_line[2] != "NA") { // place vol in its vector
			volumes.push_back(stod(tmp_line[2]));
		} else {
			volumes.push_back(nan("1"));
		}
		
		if (tmp_line[3] != "NA") { // place strike in its vector
			strikes.push_back(stod(tmp_line[3]));
		} else {
			volumes.push_back(nan("1"));
		}

		// remove the prepended & appended quotations to the option type
		string tp = tmp_line[4].substr(1, tmp_line[4].length() - 2);
		types.push_back(tp); // type

		data >> fileline; // read next line
	}
}

void OptionData::filterVolumes(double min, double max) {
	vector<string> new_dates, new_types;
	vector<double> new_lasts, new_volumes, new_strikes;
	
	for (int i = 0; i < volumes.size(); i++) {
		if ((volumes.at(i) >= min) & (volumes.at(i) <= max)) {
			new_dates.push_back( dates.at(i) );
			new_lasts.push_back( lasts.at(i) );
			new_volumes.push_back( volumes.at(i) );
			new_strikes.push_back( strikes.at(i) );
			new_types.push_back( types.at(i) );
		}
	}
	
	dates = new_dates;
	lasts = new_lasts;
	volumes = new_volumes;
	strikes = new_strikes;
	types = new_types;	
}
void OptionData::filterLasts(double min, double max) {
	vector<string> new_dates, new_types;
	vector<double> new_lasts, new_volumes, new_strikes;
	
	for (int i = 0; i < lasts.size(); i++) {
		if ((lasts.at(i) >= min) & (lasts.at(i) <= max)) {
			new_dates.push_back( dates.at(i) );
			new_lasts.push_back( lasts.at(i) );
			new_volumes.push_back( volumes.at(i) );
			new_strikes.push_back( strikes.at(i) );
			new_types.push_back( types.at(i) );
		}
	}
	
	dates = new_dates;
	lasts = new_lasts;
	volumes = new_volumes;
	strikes = new_strikes;
	types = new_types;	
}

void OptionData::computeTaus(string start_date) {
	for (int i = 0; i < dates.size(); i++) {
		taus.push_back( stringDateDiff(start_date, dates.at(i)) );
	}
}

void OptionData::computeImpliedVols(Asset A, double bisect_a, double bisect_b) {
	for (int i = 0; i < lasts.size(); i++) {
		double imp_vol;
		
		if (types.at(i) == "Call") {
			CallOption C(A, strikes.at(i), taus.at(i));
			CallOption* CO = &C;
			
			imp_vol = CO -> impliedVolBisect(lasts.at(i), bisect_a, bisect_b);
		} 
		else if (types.at(i) == "Put") {
			PutOption P(A, strikes.at(i), taus.at(i));
			PutOption* PO = &P;
						
			imp_vol = PO -> impliedVolBisect(lasts.at(i), bisect_a, bisect_b);
		} 
		else {
			cout << "WARNING in compute_implied_vol: Neither Call nor Put selected." << endl;
		}
		
		implied_vols.push_back(imp_vol);
	}
}

vector<int> OptionData::stringToDate(string str_date) { 
	// a string of format "YYYY-MM-DD" to a date
	vector<string> v_str_date;
	boost::split(v_str_date, str_date, boost::is_any_of("-"));
			
	vector<int> v_int_date(3);
	v_int_date[0] = stod(v_str_date[0]); // year
	v_int_date[1] = stod(v_str_date[1]); // month
	v_int_date[2] = stod(v_str_date[2]); // day
	
	return v_int_date;
}

tm OptionData::makeTm(vector<int> date) { 
	// make tm structure representing a date
	int year = date.at(0);
	int month = date.at(1);
	int day = date.at(2);
		
    std::tm tm = {0};
    tm.tm_year = year - 1900; // years count from 1900
    tm.tm_mon = month - 1;    // months count from January=0
    tm.tm_mday = day;         // days count from 1
    
    return tm;
}
double OptionData::stringDateDiff(string date_a, string date_b) { 
	// computes time difference in years
	vector<int> int_date_a = stringToDate(date_a);
	vector<int> int_date_b = stringToDate(date_b);
		
	tm tm_a = makeTm(int_date_a); 
	tm tm_b = makeTm(int_date_b);
		
	time_t time_a = mktime(&tm_a);
	time_t time_b = mktime(&tm_b);
	
	double tau = difftime(time_b, time_a)/(60 * 60 * 24 * 365);
		
	return tau;
}

















