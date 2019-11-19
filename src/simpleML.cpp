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

#include "simpleML.h"

using namespace std;

void simpleML::init_nodes ( node* leaf ) {
    if (leaf != 0) {
	likelihoods.insert(pair<node*,vector<double> >(leaf,static_cast<vector<double> >(0.0)));
	for (unsigned int i=1; i < n_states; ++i) likelihoods[leaf].push_back(0.0);
	init_nodes(leaf->left);
	init_nodes(leaf->right);
    }
}

void simpleML::un_init_nodes ( node* leaf ) {
    if (leaf != 0) {
       	likelihoods.erase(leaf);
        un_init_nodes(leaf->left);
        un_init_nodes(leaf->right);
    }
}

bool simpleML::check_nodes ( const node* leaf ) {
    if (leaf != 0) {
	bool return_value = true;
	if (leaf->right == 0 && leaf->left == 0) {
	    bool got_value = false;
	    for (unsigned int i=0; i < n_states; ++i)
		if (likelihoods[const_cast<node *> (leaf)][i] > 0) got_value = true;
	    if (!got_value) {
		return_value = false;
		cerr << *leaf->nodelabel << " lack character state!" << endl;
	    }
	}
	if (return_value) return_value = check_nodes(leaf->right);
        if (return_value) return_value = check_nodes(leaf->left);
	return return_value;
    }
    else return true;
}

void simpleML::calculate_likelihood (const node* leaf, sub_model& model) {
    if (leaf->left != 0) {
	calculate_likelihood (leaf->left, model);
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),leaf->left->branchlength,false, model);
    }
    if (leaf->right != 0) {
	calculate_likelihood (leaf->right, model);
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength,true, model);
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength,false,model);
    }
}

void simpleML::calculate_likelihood_rate_change_at_nodes (const node* leaf, double rate, const map<node*, double>& rate_changes, sub_model& model) {
    #ifdef DEBUG
    cerr << "Rate: " << rate << endl;
    #endif //DEBUG
    map<node*, double>::const_iterator change = rate_changes.find(const_cast<node*>(leaf));
    if (change != rate_changes.end()) rate = change->second;
    if (leaf->left != 0) {
	calculate_likelihood_rate_change_at_nodes(leaf->left, rate, rate_changes, model);
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),leaf->left->branchlength*rate,false,model);
    }
    if (leaf->right != 0) {
	calculate_likelihood_rate_change_at_nodes(leaf->right, rate, rate_changes, model);
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength*rate,true,model);
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength*rate,false,model);
    }
}

void simpleML::calculate_likelihood_rate_change_in_time(const node* leaf, const double dist_from_root, const double cut_off, double rate, sub_model& model) {
    #ifdef DEBUG
    cerr << "Rate: " << rate << "; Cut off: " << cut_off << "; Distance from root: " <<  dist_from_root << endl;
    #endif //DEBUG
    if (leaf->left != 0) {
	calculate_likelihood_rate_change_in_time(leaf->left, dist_from_root+leaf->left->branchlength,cut_off,rate, model);
	double length;
	if (dist_from_root+leaf->left->branchlength > cut_off) {
	    if (dist_from_root>cut_off)
		length = leaf->left->branchlength * rate;
	    else 
		length = (leaf->left->branchlength - (cut_off-dist_from_root))*rate + cut_off-dist_from_root;
	}
	else
	    length = leaf->left->branchlength;
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),length,false, model);
    }
    if (leaf->right != 0) {
	calculate_likelihood_rate_change_in_time(leaf->right, dist_from_root+leaf->right->branchlength,cut_off,rate, model);
	double length;
	if (dist_from_root+leaf->right->branchlength > cut_off) {
            if (dist_from_root>cut_off)
                length = leaf->right->branchlength * rate;
            else 
                length = (leaf->right->branchlength - (cut_off-dist_from_root))*rate + cut_off-dist_from_root;
        }
	else
    	    length = leaf->right->branchlength;
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),length,true,model);
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),length,false,model);
    }
}

void simpleML::branch_likelihood ( map<node*,vector<double> >::iterator LHbins, map<node*,vector<double> >::iterator startLHs, double branch_length, bool multiply, sub_model& model) {
    model.set_P_matrix(branch_length);
    for (unsigned int i = 0; i < n_states; ++i) {
            double pi = 0;
            for (unsigned int j = 0; j < n_states; ++j) {
                pi += startLHs->second[j] * model.get_P_value(i,j);
		#ifdef DEBUG
		cerr << "P " << j << " in descendant: " << startLHs->second[j] << "; P going to " << i << ": " << model.get_P_value(i,j) << endl;
		#endif //DEBUG
	    }
	    if (multiply) LHbins->second[i] *=pi;
	    else LHbins->second[i] = pi;
    }
} 

void simpleML::draw_normalized_likelihood_on_nodes( node* leaf ) {
    if ( leaf!=0 ) {
	if (leaf->left != 0 || leaf->right !=0) {
	    if (leaf->left != 0) draw_normalized_likelihood_on_nodes(leaf->left);
	    if (leaf->right !=0) draw_normalized_likelihood_on_nodes(leaf->right);
	    double sum=0;
	    for (unsigned int i=0; i < n_states; ++i) sum += likelihoods[leaf][i];
	    stringstream prob_label;
	    prob_label << "[& P={";
	    for (unsigned int i=0; i < n_states; ++i) {
		if (i!=0) prob_label << ',';
		//prob_label << "P_" << i << '=';
		prob_label << likelihoods[leaf][i]/sum;
	    }
	    prob_label << "}]";
	    string label;
	    label = prob_label.str();
	    if (leaf->nodelabel!=0) label = *leaf->nodelabel+label;
	    #ifdef DEBUG
	    std::cerr << label << endl;
	    #endif //DEBUG
	    leaf->nodelabel = nodelabels.add_string(label);
	}
    }
}

vector<character_vector> simpleML::simulate_chars( const vector<double> & charfreq, const unsigned int n_char, sub_model& model) {
    map<node*,character_vector> output_matrix;
    for (unsigned int i=0; i<n_char; ++i) { 
	map<node*,unsigned int> simdata;
	double rand_no = (double) rand()/(RAND_MAX);
	int root_trait(0);
	for (unsigned int j=0; j<charfreq.size(); ++j) {
	    rand_no-= charfreq[j];
	    if (rand_no <= 0) { root_trait = j; break; }
	}
	#ifdef DEBUG
	cerr << "Root trait: " << root_trait << endl;
	//Q_matrix.print(cerr);
	cerr << endl;
	#endif //DEBUG
	simulate_chars_subtree(root->left, simdata, root_trait, model);
	simulate_chars_subtree(root->right, simdata, root_trait, model);
	for (map<node*,unsigned int>::const_iterator k = simdata.begin(); k != simdata.end(); ++k) {
	    map<node*,character_vector>::iterator taxon = output_matrix.find(k->first);
	    bitset<SIZE> character;
	    character.set(k->second);
	    if (taxon == output_matrix.end()) {
		output_matrix[k->first] = character_vector();
		//output_matrix[k->first].set_taxon(*k->first->nodelabel);
		output_matrix[k->first].add_character(character);
	    }
	    else {
		taxon->second.add_character(character);
	    }
	}
    }
    vector<character_vector> return_data;
    for (map<node*,character_vector>::iterator i=output_matrix.begin(); i != output_matrix.end(); ++i) {
	if (i->first->nodelabel && !i->first->nodelabel->empty()) {
	    i->second.set_taxon(*i->first->nodelabel);
	    return_data.push_back(i->second);
	}
    }
    return return_data;
}

void simpleML::simulate_chars_subtree( node* leaf, map<node*,unsigned int>& simdata, const unsigned int ancestor, sub_model& model) {
    if (leaf != 0) {
	//marth::square_matrix P_matrix;
	//Q_matrix.exponential(&P_matrix, leaf->branchlength, 20);
	model.set_P_matrix(leaf->branchlength);
	unsigned int trait = ancestor;
	double rand_no = (double) rand()/(RAND_MAX);
	for (unsigned int i=0; i < model.get_n_states(); ++i) {
	    #ifdef DEBUG
	    cerr << "Rand no: " << rand_no << " P for char: " << model.get_P_value(ancestor,i) << endl; 
	    #endif //DEBUG
	    rand_no -= model.get_P_value(ancestor,i);
	    if (rand_no <= 0) { trait = i; break; }
	}
	#ifdef DEBUG
	//P_matrix.print(cerr); cerr << endl;
	if (leaf->nodelabel && !leaf->nodelabel->empty()) cerr << *leaf->nodelabel << endl;
	cerr << "Simulating trait along branch. Start: " << ancestor << " end: " << trait << endl;
	#endif //DEBUG
	simdata[leaf] = trait;
	if (leaf->left != 0) simulate_chars_subtree(leaf->left, simdata, trait, model);
	if (leaf->right != 0) simulate_chars_subtree(leaf->right, simdata, trait, model);
    }
}
