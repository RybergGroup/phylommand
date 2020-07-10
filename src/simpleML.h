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
#include "partition_models.h"
#include <map>
#include <vector>

using namespace std;

class simpleML : public tree {
    public:
	~simpleML(){
	    un_init_nodes (root);
	}
	void init() {
	    un_init_nodes (root);
	    //n_states = states;
	    init_nodes(root);
	}
	void reset() { init(); }
	void set_char( const int char_no, const string taxon, const bitset<SIZE>& state, sub_model& model ) {
	    #ifdef DEBUG
	    cerr << "Set char " << char_no << " to "<< state << " for taxon " << taxon << " (n states= " << model.get_n_states() << ")" << endl;
	    #endif //DEBUG
	    node* leaf = find_taxon_tip(root, nodelabels.find_string(taxon));
	    if (leaf != 0)
		if (char_no >= likelihoods[leaf].size()) likelihoods[leaf].resize(char_no+1);
		likelihoods[leaf][char_no].resize(model.get_n_states(),0.0);
		for (unsigned int i=0; i < model.get_n_states(); ++i) {
		    if (state[i]) likelihoods[leaf][char_no][i] = 1.0; // /state.count();
		    #ifdef DEBUG
		    cerr << "Set state " << i << " to " << likelihoods[leaf][char_no][i] << endl;
		    #endif //DEBUG
		}
	}
	double calculate_log_likelihood( sub_model& model ) {
	    #ifdef DEBUG
    	    cerr << "Calculating likelihood!!!" << endl;
	    #endif //DEBUG
	    //partitions models;
	    rate_model rate_mod;
	    return calculate_log_likelihood(model, rate_mod);
	    /*unsigned int n_char = get_n_char();
	    models.add_partition( model, rates, 0, n_char );
	    if (!check_nodes ( root )) return 0.0;
	    if (n_states != model.get_n_states()) {
		cerr << "Number of characters in model (" << model.get_n_states() << ") does not fit number of characters to explain (" << n_states << ")!!!" << endl;
		return 0.0;
	    }
	    calculate_likelihood( root, models );
	    return get_log_likelihood_from_root( models );*/
	};
	double calculate_likelihood_rate_change_at_nodes ( rate_model& rate_mod, sub_model& model ) { return calculate_log_likelihood(model, rate_mod); }
	double calculate_log_likelihood( sub_model& model, rate_model& rate_mod ) {
	    partition_models models;
	    unsigned int n_char = get_n_char( root );
	    models.add_partition( &model, &rate_mod, 0, n_char );
	    if (!check_nodes ( root, models )) return 0.0;
	    /*if (n_states != model.get_n_states()) {
                cerr << "Number of characters in model (" << model.get_n_states() << ") does not fit number of characters to explain (" << n_states << ")!!!" << endl;
                return 0.0;
            }*/
	    calculate_likelihood( root, models, 0, 0.0 );
	    #ifdef DEBUG
	    cerr << "Calculated likelihood. Now summarizing it" << endl;
	    #endif //DEBUG
            return get_log_likelihood_from_root(models);
	}
	double calculate_likelihood_rate_change_in_time(rate_model& rate_mod, sub_model& model) { return calculate_log_likelihood(model, rate_mod); }
	/*    if (!check_nodes ( root )) return 0.0;
 	    if (n_states != model.get_n_states()) {
                cerr << "Number of characters in model (" << model.get_n_states() << ") does not fit number of characters to explain (" << n_states << ")!!!" << endl;
                return 0.0;
            }
	    calculate_likelihood_rate_change_in_time( root, 0.0, rate_mod, model );
            return get_log_likelihood_from_root(model);
	}*/
	vector<character_vector> simulate_chars( const vector<double> & charfreq, const unsigned int n_char, sub_model& model);
	void draw_normalized_likelihood_on_nodes() { draw_normalized_likelihood_on_nodes( root ); }
	unsigned int get_n_char ( node* leaf );
    private:
	// variables
	map<node*,vector<vector<double> > > likelihoods;
	//vector<node*> clades;
	unsigned int n_states;
	// functions
	void simulate_chars_subtree( node* leaf, map<node*,unsigned int>& simdata, const unsigned int ancestor, sub_model& model);
	void init_nodes ( node* leaf );
	void un_init_nodes ( node* leaf );
	bool check_nodes ( const node* leaf, partition_models& models );
	void calculate_likelihood(const node* leaf, partition_models& partitions, unsigned int clade_no, const double dist_from_root);
	/*void calculate_likelihood(const node* leaf, sub_model& model, rate_model& rate_mod, unsigned int clade_no, const double dist_from_root);
	void calculate_likelihood ( const node* leaf, sub_model& model ) {
	    rate_model rate_mod;
	    calculate_likelihood(leaf, model, rate_mod, 0, 0);
	}
	void calculate_likelihood_rate_change_in_time(const node* leaf, const double dist_from_root, rate_model& rate_mod, sub_model& model) { calculate_likelihood(leaf, model, rate_mod, 0, dist_from_root); }
	void calculate_likelihood_rate_change_at_nodes(const node* leaf, unsigned int clade_no, rate_model& rate_mod, sub_model& model) { calculate_likelihood(leaf, model, rate_mod, clade_no, 0 ); }*/
	void draw_normalized_likelihood_on_nodes( node* leaf );
	void branch_likelihood ( vector<double>& LHbins, vector<double>& startLHs, double branch_length, bool multiply, sub_model& model);
	double get_log_likelihood_from_root ( partition_models& models ) {
	    double log_likelihood(0.0);
	    #ifdef DEBUG
	    cerr << "Starting to summerize for " << likelihoods[root].size() << " positions" << endl;
	    #endif //DEBUG
	    for (unsigned int pos = 0; pos < likelihoods[root].size(); ++pos) {
		sub_model* model  = models.get_sub_model( pos );
		double likelihood(0.0);
		for ( unsigned int char_no=0; char_no < model->get_n_states() && char_no < likelihoods[root][pos].size(); ++char_no ) {
		    if ( model->is_tr() ) likelihood += likelihoods[root][pos][char_no] * model->get_freq(char_no);
		    else likelihood += likelihoods[root][pos][char_no];
		}
		log_likelihood += log(likelihood);
		#ifdef DEBUG
		cerr << "Log like char " << pos << " is " << log_likelihood << endl;
		#endif //DEBUG
	    }
	    return log_likelihood;
	}
};
