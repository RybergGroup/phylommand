#include "sub_model.h"

unsigned int sub_model::get_rate_pos (const unsigned int row, const unsigned int col) { // returns the number of possible rates if not finding possition, i.e. out of bound
    if (row==col) return n_posible_rates();
    unsigned int pos(0);
    if (timerev) {
	unsigned int i(0);
	unsigned int j(0);
	if (col>row) { // if in right hand of matrix
	    i=row; // go by row
	    j=col-row; // column minus columns in the left of the matrix
	}
	else if (row>col) { // same but by column
	    i=col;
	    j=row-col;
	}
	while(i) { pos += n_states-i; --i; }
	pos += j-1;
    }
    else {
	for (unsigned int i=0; i < row; ++i) pos += n_states-1; // add n_states - diagonal for each row
	pos += col; // add the number of columns
	if (col > row) --pos; // but if passing diagonal remove it
    }
    return pos;
}

bool sub_model::set_freq (const unsigned int n, const double value) {
    if (n >= n_states || value < 0.0) return false;
    double sum(0.0);
    for (unsigned int i=0; i < n_states-1; ++i) {
	if (i==n) sum += value;
	else sum += state_freqs[i];
    }
    if (sum > 1.0) return false;
    state_freqs[n] = value;
    calc_Q_matrix();
    return true;
}

double sub_model::get_freq (const unsigned int n) {
    if (n >= n_states) return 0.0;
    else if (equal_freq) return 1.0/n_states;
    else if (n == n_states-1) {
	double freq(1.0);
	for (unsigned int i=0; i < n_states-1; ++i) freq -= state_freqs[i];
	return freq;
    }
    else return state_freqs[n];
}

bool sub_model::set_rate_spec (const unsigned int n, const unsigned int param) {
    // test that n is within bounds
    const unsigned int n_pos_rates(n_posible_rates());
    if (n > n_pos_rates) return false; // if higher than number of possible parameters
    // test so all lower parameter values are set
    bitset<400> numbers;
    for (unsigned int i = 0; i < n_pos_rates; ++i) {
	if (i == n) continue;
	numbers.set(rate_specs[i]);
    }
    unsigned int max(0);
    for (unsigned int i = 0; i < n_pos_rates; ++i) {
	if (i < param && !numbers.test(i)) return false;
	if (numbers.test(i) && i > max) max = i;
    }
    rate_specs[n] = param;
    // Make room for new parameter if needed
    if (max+1 > n_rates) {
	double* temp = new double(max+1);
	for (unsigned int i=0; i < n_rates; ++i) temp[i] = rates[i];
	delete rates;
	rates = temp;
	n_rates = max+1;
    }
    return true;
}

void sub_model::calc_Q_matrix() {
    for (unsigned int i=0; i<n_states; ++i) {
	double neg_sum(0.0);
	for (unsigned int j=0; j < n_states; ++j) {
	    if (i==j) continue;
	    double value = get_rate(i,j);
	    if (timerev) value *= get_freq(j);
	    neg_sum -= value;
	    Q_matrix->set_value(i,j,value);
	}
	Q_matrix->set_value(i,i,neg_sum);
    }
}
void sub_model::set_P_matrix( const double branch_length ) {
    if (P_matrix == 0) P_matrix = new marth::square_matrix(n_states);
    else {
	delete P_matrix;
	P_matrix = new marth::square_matrix(n_states);
    }
    Q_matrix->exponential(P_matrix,branch_length,20);
}
