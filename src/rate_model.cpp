/********************************************************************
Copyright (C) 2020 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

contact: martin.ryberg@ebc.uu.se
*********************************************************************/

#include "rate_model.h"

using namespace std;

bool rate_model::add_rate_for_time (unsigned int rate_class, const float time) {
    if (rate_class >= n_rates_in_time) return false;
    if (rate_time_cut_offs == 0) {
	rate_time_cut_offs = new rate_cut_off[1];
	rate_time_cut_offs[0].cut_off = time;
	rate_time_cut_offs[0].rate_class = rate_class;
	++n_rate_time_cut_offs;
	return true;
    }
    rate_cut_off* new_cut_offs = new rate_cut_off[n_rate_time_cut_offs+1];
    unsigned int new_pos(0);
    unsigned int i(0);
    for (; i < n_rate_time_cut_offs; ++i) {
        if (new_pos == i && rate_time_cut_offs[i].cut_off > time) {
	    new_cut_offs[new_pos].cut_off = time;
	    new_cut_offs[new_pos].rate_class = rate_class;
	    ++new_pos;
	}
	new_cut_offs[new_pos].cut_off = rate_time_cut_offs[i].cut_off;
	new_cut_offs[new_pos].rate_class = rate_time_cut_offs[i].rate_class;
	++new_pos;
    }
    if (new_pos == i) {
	new_cut_offs[new_pos].cut_off = time;
	new_cut_offs[new_pos].rate_class = rate_class;
    }
    delete[] rate_time_cut_offs;
    rate_time_cut_offs = new_cut_offs;
    ++n_rate_time_cut_offs;
    return true;
}

void rate_model::set_time_rate (const unsigned int rate_class, const float rate) {
    if (rate_class >= n_rates_in_time) {
	float* new_rates = new float[rate_class+1];
	unsigned int i(0);
	for (i = 0; i < n_rates_in_time; ++i) new_rates[i] = rates_in_time[i];
	for (; i <=rate_class; ++i) new_rates[i] = rate;
	if (rates_in_time != 0) delete[] rates_in_time;
	rates_in_time = new_rates;
	n_rates_in_time = rate_class+1;
    }
    else rates_in_time[rate_class] = rate;
}
bool rate_model::set_clade_rate ( const unsigned int clade, const unsigned int rate_class ) {
    // add checks
    if (rate_class >= n_rates) return false;
    else if (clade >= n_clade_rates) {
	unsigned int* new_clade_rates = new unsigned int[clade+1];
	unsigned int i(0);
	for (; i < n_clade_rates; ++i) new_clade_rates[i] = clade_rates[i];
	for (; i <= clade; ++i) new_clade_rates[i] = rate_class;
	if (clade_rates != 0) delete[] clade_rates;
	clade_rates = new_clade_rates;
	n_clade_rates = clade+1;
    }
    else clade_rates[clade] = rate_class;
    return true;
};
void rate_model::set_rate (const unsigned int rate_class, const float rate) {
    if (rate_class >= n_rates) {
	float* new_rates = new float[rate_class+1];
	unsigned int i(0);
	for (; i < n_rates; ++i) new_rates[i] = rates[i];
	for (; i <= rate_class; ++i) new_rates[i] = rate;
	if (rates != 0) delete[] rates;
	rates = new_rates;
	n_rates = rate_class+1;
    }
    else rates[rate_class] = rate;
}

bool rate_model::set_parameter (const unsigned int n, const float rate) {
    if (n < n_rate_time_cut_offs) return set_cut_off(n,rate);
    else if (n-n_rate_time_cut_offs < n_rates_in_time) { set_time_rate(n-n_rate_time_cut_offs,rate); return true; }
    else if (n-n_rate_time_cut_offs-n_rates_in_time < n_rates) { set_rate(n-n_rate_time_cut_offs-n_rates_in_time, rate); return true; }
    return false;
}

double rate_model::get_time_mod_branch_length (const double branch_length, const float distance_to_root) {
    const double branch_end = distance_to_root+branch_length;
    double new_br_length(0.0);
    double previous_cut_off(-0.0);
    for (unsigned int i=0; i < n_rate_time_cut_offs; ++i) {
        if ( distance_to_root < rate_time_cut_offs[i].cut_off && branch_end > previous_cut_off) { //
            if (distance_to_root >= previous_cut_off && branch_end <= rate_time_cut_offs[i].cut_off) // if entirely within interval
                new_br_length += (branch_end-distance_to_root)*rates_in_time[rate_time_cut_offs[i].rate_class]; // multiply entire branch
            else if (distance_to_root < previous_cut_off && branch_end <= rate_time_cut_offs[i].cut_off) // if also within previous interval
                new_br_length += (branch_end-previous_cut_off)*rates_in_time[rate_time_cut_offs[i].rate_class]; // multiply what is left of branch
            else if (distance_to_root >= previous_cut_off && branch_end > rate_time_cut_offs[i].cut_off) // if also in next interval
                new_br_length += (rate_time_cut_offs[i].cut_off-distance_to_root)*rates_in_time[rate_time_cut_offs[i].rate_class]; // multiply from start of branch until end of interval
            else if (distance_to_root < previous_cut_off && branch_end > rate_time_cut_offs[i].cut_off) // if in also in both previous and next interval
                new_br_length += (rate_time_cut_offs[i].cut_off-previous_cut_off)*rates_in_time[rate_time_cut_offs[i].rate_class]; // multiply the distance between the intervals
        }
        previous_cut_off = rate_time_cut_offs[i].cut_off;
    }
    if (distance_to_root < previous_cut_off && branch_end > previous_cut_off) // if sticking out beyond last interval
        new_br_length += branch_end-previous_cut_off; // add what is left
    else if (distance_to_root > previous_cut_off) new_br_length = branch_length;
/*    #ifdef DEBUG
    cerr << "Branch start: " << distance_to_root << " Branch end: " << branch_end << " Branch length in: " << branch_length << " Branch length out: " << new_br_length << endl;
    #endif //DEBUG */
    return new_br_length;
}

bool rate_model::pars_clade_rates_arg(const char* argument, vector<string>& taxon_vector) {
    unsigned int rate_class(0);
    string temp;
    for (unsigned int j=0; argument[j] != '\0'; ++j) {
        if (argument[j] == ':') {
            rate_class = get_n_clade_rates();
            set_rate(rate_class,atof(temp.c_str()));
            temp.clear();
        }
        else if (argument[j] == ';' || argument[j+1] == '\0') {
	    if (argument[j+1] == '\0') temp += argument[j];
            taxon_vector.push_back(temp);
	    #ifdef DEBUG
	    cerr << "Seting rate (" << rate_class << "): " << rates[rate_class] << " for clade " << taxon_vector.size()-1 << " \"" << taxon_vector[taxon_vector.size()-1] << "\"" << endl;
	    #endif //DEBUG
            if (!set_clade_rate(taxon_vector.size()-1,rate_class)) {
                set_rate(rate_class,1.0);
                if (!set_clade_rate(taxon_vector.size()-1,rate_class)) return false;
            }
            temp.clear();
        }
        else temp += argument[j];
    }
    if (get_n_clade_rates() > 0) return true;
    else return false;
}

bool rate_model::pars_rates_and_times_arg(const char* argument) {
    string temp;
    //argv_parser::pars_sub_args(argv[i], ',', arguments );
    unsigned int rate_class(0);
    for (unsigned int i=0; argument[i] != '\0'; ++i) {
	if (argument[i] == ':') {
	    rate_class = get_n_time_rates();
	    set_time_rate(rate_class, atof(temp.c_str()));
	    temp.clear();
	}
	else if (argument[i] == ',' || argument[i+1] == '\0') {
	    if (argument[i+1] == '\0') temp += argument[i];
	    if (!add_rate_for_time(rate_class,atof(temp.c_str()))) return false;
	    temp.clear();
	}
	else temp += argument[i];
    }
    if (get_n_time_rates() > 0) return true;
    else return false;
}

