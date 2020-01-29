/********************************************************************
Copyright (C) 2020 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

contact: martin.ryberg@ebc.uu.se
*********************************************************************/

#include <iostream>
#include "tree.h"
#include <cmath>
#include "sub_model.h"
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
	    init_nodes(root);
	};
	void reset(const unsigned int states) { init(states); };
	void set_char( const string taxon, const bitset<SIZE>& state) {
	    node* leaf = find_taxon_tip(root, nodelabels.find_string(taxon));
	    if (leaf != 0)
		for (unsigned int i=0; i < n_states; ++i) {
		    if (state[i]) likelihoods[leaf][i] = 1.0; // /state.count();
		    else likelihoods[leaf][i] = 0.0;
		}
	};
	double calculate_log_likelihood( sub_model& model ) {
	    #ifdef DEBUG
    	    cerr << "Calculating likelihood!!!" << endl;
	    #endif //DEBUG
	    if (!check_nodes ( root )) return 0.0;
	    if (n_states != model.get_n_states()) {
		cerr << "Number of characters in model (" << model.get_n_states() << ") does not fit number of characters to explain (" << n_states << ")!!!" << endl;
		return 0.0;
	    }
	    calculate_likelihood( root, model );
	    double likelihood=0;
	    for (unsigned int i=0; i < n_states; ++i) likelihood += likelihoods[root][i];
	    return log(likelihood);
	};
	double calculate_likelihood_rate_change_at_nodes ( rate_model& rate_mod, sub_model& model ) {
	    if (!check_nodes ( root )) return 0.0;
	    if (n_states != model.get_n_states()) {
                cerr << "Number of characters in model (" << model.get_n_states() << ") does not fit number of characters to explain (" << n_states << ")!!!" << endl;
                return 0.0;
            }
	    calculate_likelihood_rate_change_at_nodes(root,1.0,rate_mod, model);
	    double likelihood=0;
	    for (unsigned int i=0; i < n_states; ++i) likelihood += likelihoods[root][i];
            return log(likelihood);
	}
	double calculate_likelihood_rate_change_in_time(rate_model& rate_mod, sub_model& model) {
	    if (!check_nodes ( root )) return 0.0;
 	    if (n_states != model.get_n_states()) {
                cerr << "Number of characters in model (" << model.get_n_states() << ") does not fit number of characters to explain (" << n_states << ")!!!" << endl;
                return 0.0;
            }
	    calculate_likelihood_rate_change_in_time( root, 0.0, rate_mod, model );
	    double likelihood=0;
            for (unsigned int i=0; i < n_states; ++i) likelihood += likelihoods[root][i];
            return log(likelihood);
	}
	vector<character_vector> simulate_chars( const vector<double> & charfreq, const unsigned int n_char, sub_model& model);
	void draw_normalized_likelihood_on_nodes() { draw_normalized_likelihood_on_nodes( root ); };
	unsigned int n_taxon_sets() { return clades.size(); }
    private:
	// variables
	map<node*,vector<double> > likelihoods;
	//vector<node*> clades;
	unsigned int n_states;
	// functions
	void simulate_chars_subtree( node* leaf, map<node*,unsigned int>& simdata, const unsigned int ancestor, sub_model& model);
	void init_nodes ( node* leaf );
	void un_init_nodes ( node* leaf );
	bool check_nodes ( const node* leaf );
	void calculate_likelihood ( const node* leaf, sub_model& model );
	void calculate_likelihood_rate_change_in_time(const node* leaf, const double dist_from_root, rate_model& rate_mod, sub_model& model);
	void calculate_likelihood_rate_change_at_nodes(const node* leaf, double rate, rate_model& rate_mod, sub_model& model);
	void draw_normalized_likelihood_on_nodes( node* leaf );
	void branch_likelihood ( map<node*,vector<double> >::iterator LHbins, map<node*,vector<double> >::iterator startLHs, double branch_length, bool multiply, sub_model& model);
};
