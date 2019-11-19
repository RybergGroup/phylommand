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
#ifndef SUB_MODELHEADER
#define SUB_MODELHEADER

#include "marth/marth.h"
#include <vector>
#include <bitset>

using namespace std;

class sub_model {
    public:
    sub_model ( unsigned int dim, bool tr, bool eq ) : timerev(tr), equal_freq(eq), n_states(dim), n_rates(1), state_freqs(0), P_matrix(0) {
	set_sub_model(eq);
    }
    sub_model ( unsigned int dim, bool tr = false) : timerev(tr), n_states(dim), n_rates(1), state_freqs(0), P_matrix(0) {
	set_sub_model( false );
	#ifdef DEBUG
	cerr << "Number of dimentions: " << dim << " number of states: " << n_states << endl;
	#endif //DEBUG
    }
    ~sub_model() {
	if (rate_specs != 0) delete [] rate_specs;
	if (rates != 0) delete [] rates;
	if (state_freqs != 0) delete [] state_freqs;
	if (P_matrix != 0) delete P_matrix;
	delete Q_matrix;
    }
    void set_P_matrix( const double branch_length );
    double get_P_value( const unsigned int row, const unsigned int col) {
	if (row < n_states && col < n_states) return P_matrix->get_value(row,col);
	else return 0.0;
    }
    bool set_rate (const unsigned int n, const double value) {
	if (n < n_rates && value >= 0.0) {
	    rates[n] = value;
	    calc_Q_matrix();
	    return true;
	}
	else return false;
    }
    bool set_rate (const unsigned int row, const unsigned int col, const double value) {
	unsigned int pos = get_rate_pos( row, col );
	if (pos < n_posible_rates())
	    return set_rate(rate_specs[pos],value);
	else return false;
    }
    double get_rate(const unsigned int n) {
	if (n < n_rates) return rates[n];
	else return 0.0;
    }
    double get_rate(const unsigned int row, const unsigned int col) {
	unsigned int pos = get_rate_pos( row, col );
	if (pos < n_posible_rates())
	    return rates[rate_specs[pos]];
	else return 0.0;
    }
    bool set_rate_spec (const unsigned int n, const unsigned int param);
    bool set_rate_spec (const unsigned int row, const unsigned int col, const unsigned int param) {
	unsigned int pos = get_rate_pos( row, col );
	return set_rate_spec ( pos, param );
    }
    unsigned int get_n_rates () { return n_rates; }
    bool set_freq (const unsigned int n, const double value);
    void set_freq_same () {
	if (state_freqs != 0)
	    for (unsigned int i=0; i < n_states-1; ++i) state_freqs[i] = 1.0/n_states;
    }
    double get_freq (const unsigned int n);
    void equal_state_freq ( const bool eq) { equal_freq=eq; }
    bool is_time_reversible() { return timerev; }
    unsigned int get_n_states () { return n_states; }
    bool is_timerev() { return timerev; }
    bool is_tr() { return timerev; }
    bool are_state_freq_eqal() { return equal_freq; }
    unsigned int n_posible_rates() {
	if (timerev) return ((n_states*n_states)-n_states)/2;
	else return (n_states*n_states)-n_states;
    }
    void print_Q_matrix ( ostream& output_stream ) { Q_matrix->print( output_stream ); }
    void print_Q_matrix () { Q_matrix->print( std::cout ); }
    private:
    void set_sub_model( bool eq ) {
	equal_freq = eq;
	if (timerev) {
	    rate_specs = new unsigned int[((n_states*n_states)-n_states)/2];
	    for (unsigned int i=0; i < ((n_states*n_states)-n_states)/2; ++i) rate_specs[i] = 0;
	}
	else {
	    rate_specs = new unsigned int[(n_states*n_states)-n_states];
	    for (unsigned int i=0; i < ((n_states*n_states)-n_states); ++i) rate_specs[i] = 0;
	}
	rates = new double[1];
	Q_matrix = new marth::square_matrix(n_states);
    }
    const bool timerev;
    bool equal_freq;
    const unsigned int n_states;
    unsigned int n_rates;
    unsigned int* rate_specs;
    double* rates;
    double* state_freqs;
    marth::square_matrix* Q_matrix;
    marth::square_matrix* P_matrix;
    unsigned int get_rate_pos (const unsigned int row, const unsigned int col);
    void calc_Q_matrix();
};

#endif //SUB_MODELHEADER
