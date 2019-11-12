/********************************************************************
Copyright (C) 2019 Martin Ryberg

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

#include "marth/marth.h"
#include <vector>
#include <bitset>

using namespace std;

class sub_model {
    public:
    sub_model( unsigned int dim, bool tr, bool eq ) : timerev(tr), equal_freq(eq), n_states(dim), state_freqs(0), n_rates(1) {
	if (timerev) {
	    rates_specs = new unsigned int[((n_states*n_states)-n_states)/2];
	}
	else {
	    rates_specs = new unsigned int[(n_states*n_states)-n_states]
	}
	rates = new double[1];
	Q_matrix = new marth::square_matrix(n_states);
    }
    sub_model( unsigned int dim, bool tr) { sub_model(dim, tr, false); }
    sub_model( unsigned int dim ) { sub_model(dim, false, false); }
    ~sub_model() {
	if (rates_specs != 0) delete [] rates_specs;
	if (freq_specs != 0) delete [] freq_specs;
	if (rates != 0) delete [] rates;
	if (state_freqs != 0) delete [] state_freqs;
	delete Q_matrix;
    }
    marth::square_matrix get_P_matrix( double branch );
    void set_rate (const unsigned int n, const double value) { if (n < n_rates) rates[n] = value; }
    void set_rate (const unsigned int row, const unsigned int col, const double value) {
	unsigned int pos = get_rate_pos( row, col );
	if (pos < n_posible_rates())
	    set_rate(rate_spec[pos],value);
    }
    bool set_rate_spec (const unsigned int n, const unsigned int param);
    bool set_rate_spec (const unsigned int row, const unsigned int col, const unsigned int param) {
	unsigned int pos = get_rate_pos( row, col );
	return set_rate_spec ( pos, param );
    }
    void set_freq (unsigned int n, double value);
    double get_freq (unsigned int n);
    void equal_state_freq ( bool eq) { equal_freq=eq; }
    bool is_time_reversible() { return timerev; }
    bool is_timerev() { return timerev; }
    bool is_tr() () { return timerev; }
    unsigned int n_distinct_freqs ();
    unsigned int n_posible_rates() {
	if (timerev) return ((n_states*n_states)-n_states)/2;
	else return (n_states*n_states)-n_states;
    }
    private:
    bool timerev;
    bool equal_freq;
    unsigned int n_states;
    unsigned int n_rates;
    unsigned int* rate_specs;
    double* rates;
    double* state_freqs;
    marth::square_matrix* Q_matrix;
    get_rate_pos (unsigned int row, unsigned int col);

}
