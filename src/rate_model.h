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

#ifndef RATEMODHEADER
#define RATEMODHEADER

#include <iostream>
#include <string.h>
#include <vector>

using namespace std;

class rate_model {
    public:
	rate_model(): rate_time_cut_offs(0), n_rate_time_cut_offs(0), rates_in_time(0), n_rates_in_time(0), rates(0), n_rates(0), clade_rates(0), n_clade_rates(0) { };
	~rate_model() {
	    if (rate_time_cut_offs != 0) delete[] rate_time_cut_offs;
	    if (rates_in_time != 0) delete[] rates_in_time;
	    if (rates != 0) delete[] rates;
	    if (clade_rates != 0) delete[] clade_rates;
	}
    	bool add_rate_for_time (const unsigned int rate_class, const float time);
	bool set_cut_off (const unsigned int cut_off, float value) {
	    if (cut_off > n_rate_time_cut_offs) return false;
	    if (cut_off > 0 && value < rate_time_cut_offs[cut_off-1].cut_off) return false;
	    if (cut_off < n_rate_time_cut_offs-1 && value > rate_time_cut_offs[cut_off+1].cut_off) return false;
	    rate_time_cut_offs[cut_off].cut_off = value;
	    return true;
	}
        void set_time_rate (const unsigned int rate_class, const float rate);
        bool set_clade_rate ( const unsigned int clade, const unsigned int rate_class );
        void set_rate (const unsigned int rate_class, const float rate);
	bool set_parameter (const unsigned int n, const float rate);
        unsigned int get_n_rates () { return n_rates_in_time+n_rates; }
        unsigned int get_n_clade_rates () { return n_rates; }
        unsigned int get_n_time_rates () { return n_rates_in_time; }
	unsigned int get_n_cut_offs () { return n_rate_time_cut_offs; }
	unsigned int get_n_clades () { return n_clade_rates; }
	unsigned int get_n_parameters() { return n_rates_in_time+n_rates+n_rate_time_cut_offs; }
	float get_cut_off (const unsigned int n) {
	    if (n >= n_rate_time_cut_offs) return 0.0;
	    return rate_time_cut_offs[n].cut_off;
	}
	float get_rate_in_time (const unsigned int n) {
	    if (n >= n_rate_time_cut_offs) return 1.0;
	    return rates_in_time[rate_time_cut_offs[n].rate_class];
	}
        float get_rate_clade( const unsigned int n ) {
            if (n < n_clade_rates) return rates[clade_rates[n]];
            else return 1;
        }
        double get_time_mod_branch_length (const double branch_length, const float distance_to_root);
	bool pars_clade_rates_arg (const char* argument, vector<string>& taxon_string);
	bool pars_rates_and_times_arg (const char* argument);
    private:
	struct rate_cut_off {
	    unsigned int rate_class;
	    float cut_off;
	};
	rate_cut_off* rate_time_cut_offs;
	unsigned int n_rate_time_cut_offs;
        float* rates_in_time;
	unsigned int n_rates_in_time;
        float* rates;
	unsigned int n_rates;
        unsigned int* clade_rates;
	unsigned int n_clade_rates;
};

#endif //RATEMODHEADER
