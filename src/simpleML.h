/********************************************************************
Copyright (C) 2016 Martin Ryberg

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

#include <iostream>
#include "tree.h"
#include <cmath>
#include "marth/marth.h"
#include <map>
#include <vector>

using namespace std;

class simpleML : public tree {
    public:
	~simpleML(){
	    un_init_nodes (root);
	}
	void init(unsigned int states) {
	    un_init_nodes (root);
	    n_states = states;
	    Q_matrix.reset(n_states);
	    init_nodes(root);
	};
	void reset(const unsigned int states) { init(states); };
	void set_char( const string taxon, const bitset<SIZE>& state) { //const double* state ) {
	    node* leaf = find_taxon_tip(root, nodelabels.find_string(taxon));
	    if (leaf != 0)
		for (unsigned int i=0; i < n_states; ++i) {
		    if (state[i]) likelihoods[leaf][i] = 1.0; //state[i];
		    else likelihoods[leaf][i] = 0.0;
		}
	};
	double calculate_log_likelihood() {
	    if (!check_nodes ( root )) return 0.0;
	    calculate_likelihood( root );
	    double likelihood=0;
	    for (unsigned int i=0; i < n_states; ++i) likelihood += likelihoods[root][i];
	    return log(likelihood);
	};
	void set_Q_matrix ( const unsigned int* parameters, const double* values );
	void set_Q_matrix ( const double* values );
	double calculate_log_likelihood( const unsigned int* parameters, const double* values ) {
	    set_Q_matrix(parameters,values);
	    return calculate_log_likelihood();
	};
	double calculate_log_likelihood( const double* values ) {
	    set_Q_matrix(values);
	    return calculate_log_likelihood();
	}
	void add_taxon_sets(vector<string>& taxon_sets) {
	    for (int i = 0; i < taxon_sets.size(); ++i) {
		clades.push_back(most_recent_common_ancestor(taxon_sets[i]));
	    }
	}
	double calculate_likelihood_rate_change_at_nodes ( const double* rates ) {
	    map<node*, double> rate_changes;
	    for (int i = 0; i < clades.size() && i < clades.size(); ++i) {
		if (clades[i] != 0) rate_changes[clades[i]] = rates[i];
	    }
	    calculate_likelihood_rate_change_at_nodes(root,1.0,rate_changes);
	    double likelihood=0;
	    for (unsigned int i=0; i < n_states; ++i) likelihood += likelihoods[root][i];
            return log(likelihood);
	}
	double calculate_likelihood_rate_change_at_nodes (const double* values, const double* rates) {
	    set_Q_matrix(values);
	    return calculate_likelihood_rate_change_at_nodes( rates);
	}
	double calculate_likelihood_rate_change_at_nodes (const unsigned int* parameters, const double* values, const double* rates) {
	    set_Q_matrix(parameters,values);
	    return calculate_likelihood_rate_change_at_nodes( rates);
	}
	double calculate_likelihood_rate_change_in_time(const double cut_off, const double rate) {
	    if (!check_nodes ( root )) return 0.0;
	    calculate_likelihood_rate_change_in_time( root, 0.0, cut_off, rate );
	    double likelihood=0;
            for (unsigned int i=0; i < n_states; ++i) likelihood += likelihoods[root][i];
            return log(likelihood);
	}
	double calculate_likelihood_rate_change_in_time(const unsigned int* parameters, const double* values, const double cut_off, const double rate ) {
	    set_Q_matrix(parameters,values);
	    return calculate_likelihood_rate_change_in_time(cut_off, rate);
	}
	double calculate_likelihood_rate_change_in_time( const double* values, const double cut_off, const double rate ) {
	    set_Q_matrix(values);
	    return calculate_likelihood_rate_change_in_time(cut_off, rate);
        }
	void print_Q_matrix ( ostream& output_stream ) { Q_matrix.print( output_stream ); };
	void print_Q_matrix () { Q_matrix.print( std::cout ); };
	void draw_normalized_likelihood_on_nodes() { draw_normalized_likelihood_on_nodes( root ); };
	unsigned int n_taxon_sets() { return clades.size(); }
    private:
	// variables
	map<node*,vector<double> > likelihoods;
	vector<node*> clades;
	unsigned int n_states;
	marth::square_matrix Q_matrix;
	// functions
	void init_nodes ( node* leaf );
	void un_init_nodes ( node* leaf );
	bool check_nodes ( const node* leaf );
	void calculate_likelihood ( const node* leaf );
	void calculate_likelihood_rate_change_in_time(const node* leaf, const double dist_from_root, const double cut_off, const double rate);
	void calculate_likelihood_rate_change_at_nodes(const node* leaf, double rate, const map<node*, double>& rate_changes);
	void draw_normalized_likelihood_on_nodes( node* leaf );
	void branch_likelihood ( map<node*,vector<double> >::iterator LHbins, map<node*,vector<double> >::iterator startLHs, double branch_length, bool multiply);
};
