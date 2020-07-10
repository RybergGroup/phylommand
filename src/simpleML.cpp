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
	likelihoods.insert(pair<node*,vector<vector<double> > >(leaf,static_cast<vector<vector<double> > >(0.0)));
	//for (unsigned int i=1; i < n_states; ++i) likelihoods[leaf].push_back(0.0);
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

unsigned int simpleML::get_n_char ( node* leaf ) {
    if (leaf == 0) return 0;
    unsigned int left(0);
    unsigned int right(0);
    if (leaf->left == 0 && leaf->right ==0 && likelihoods.find(leaf) != likelihoods.end()) return likelihoods[const_cast<node *> (leaf)].size(); 
    if (leaf->left != 0) left = get_n_char(leaf->left);
    if (leaf->right != 0) right = get_n_char(leaf->right);
    if (left > right) return left;
    return right;
}

bool simpleML::check_nodes ( const node* leaf, partition_models& models ) {
    if (leaf != 0) {
	bool return_value = true;
	if (leaf->right == 0 && leaf->left == 0) {
	    for (unsigned int pos=0; pos < likelihoods[const_cast<node *> (leaf)].size(); ++pos) {
		bool got_value = false;
		for (unsigned int i=0; i < models.get_n_states(pos); ++i)
		    if (likelihoods[const_cast<node *> (leaf)][pos][i] > 0) got_value = true;
		if (!got_value) {
		    return_value = false;
		    cerr << *leaf->nodelabel << " lack character state for character" << pos << "!" << endl;
		}
	    }
	}
	if (return_value) return_value = check_nodes(leaf->right, models);
        if (return_value) return_value = check_nodes(leaf->left, models);
	return return_value;
    }
    else return true;
}

/*void simpleML::calculate_likelihood (const node* leaf, sub_model& model) {
    if (leaf->left != 0) {
	calculate_likelihood (leaf->left, model);
	map<node*,vector<vector<double> > >::iterator node_lik = likelihoods.find(const_cast<node*>(leaf));
	map<node*,vector<vector<double> > >::iterator node_lik_left = likelihoods.find(const_cast<node*>(leaf->left));
	if (node_lik != likelihoods.end() && node_lik_left != likelihoods.end()) {
	    for (unsigned int pos = 0; pos < node_lik.size(); ++pos) {
		branch_likelihood(node_lik->at(pos),node_lik_left->at(pos),leaf->left->branchlength,false, model.get_model(pos));
	    }
	}
    }
    if (leaf->right != 0) {
	calculate_likelihood (leaf->right, model);
	map<node*,vector<vector<double> > >::iterator node_lik = likelihoods.find(const_cast<node*>(leaf));
	map<node*,vector<vector<double> > >::iterator node_lik_right = likelihoods.find(const_cast<node*>(leaf->right));
	if (node_lik != likelihoods.end() && node_lik_left != likelihoods.end()) {
	    for (unsigned int pos = 0; pos < node_lik.size(); ++pos) {
		branch_likelihood(node_lik->at(pos),node_lik_left->at(pos),leaf->left->branchlength,false, model.get_model(pos));
		if (leaf->left != 0) branch_likelihood(node_lik->at(pos),node_lik_right->at(pos),leaf->right->branchlength,true, model.get_model(0));
		else branch_likelihood(likelihoods.find(node_lik->at(pos),node_lik_left->at(pos),leaf->right->branchlength,false,model.get_model(0));
	    }
	}
    }
}

void simpleML::calculate_likelihood_rate_change_at_nodes (const node* leaf, double rate, rate_model& rate_mod, sub_model& model) {
    #ifdef DEBUG
    cerr << "Rate: " << rate << endl;
    #endif //DEBUG
    map<node*, unsigned int>::const_iterator change = clades.find(const_cast<node*>(leaf));
    if (change != clades.end()) rate = rate_mod.get_rate_clade(change->second);
    if (leaf->left != 0) {
	calculate_likelihood_rate_change_at_nodes(leaf->left, rate, rate_mod, model);
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),leaf->left->branchlength*rate,false,model.get_model(0));
    }
    if (leaf->right != 0) {
	calculate_likelihood_rate_change_at_nodes(leaf->right, rate, rate_mod, model);
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength*rate,true,model.get_model(0));
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength*rate,false,model.get_model(0));
    }
}

void simpleML::calculate_likelihood_rate_change_in_time(const node* leaf, const double dist_from_root, rate_model& rate_mod, sub_model& model) {
    if (leaf->left != 0) {
	calculate_likelihood_rate_change_in_time(leaf->left, dist_from_root+leaf->left->branchlength, rate_mod, model);
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),rate_mod.get_time_mod_branch_length(leaf->left->branchlength,dist_from_root),false, model.get_model(0));
    }
    if (leaf->right != 0) {
	calculate_likelihood_rate_change_in_time(leaf->right, dist_from_root+leaf->right->branchlength,rate_mod, model);
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),rate_mod.get_time_mod_branch_length(leaf->right->branchlength,dist_from_root),true,model.get_model(0));
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),rate_mod.get_time_mod_branch_length(leaf->right->branchlength,dist_from_root),false,model.get_model(0));
    }
}*/

void simpleML::calculate_likelihood(const node* leaf, partition_models& partitions, unsigned int clade_no, const double dist_from_root) {
    if (leaf->left == 0 && leaf->right == 0) return;
    map<node*, unsigned int>::const_iterator change = clades.find(const_cast<node*>(leaf));
    map<node*,vector<vector<double> > >::iterator node_lik_left;
    map<node*,vector<vector<double> > >::iterator node_lik_right;
    if (change != clades.end()) clade_no = change->second;
    if (leaf->left != 0) {
	calculate_likelihood(leaf->left,partitions,clade_no,dist_from_root+leaf->left->branchlength);
	node_lik_left = likelihoods.find(const_cast<node*>(leaf->left));
    }
    if (leaf->right != 0) {
	calculate_likelihood(leaf->right,partitions,clade_no,dist_from_root+leaf->right->branchlength);
	node_lik_right = likelihoods.find(const_cast<node*>(leaf->right));
    }
    map<node*,vector<vector<double> > >::iterator node_lik = likelihoods.find(const_cast<node*>(leaf));
    if (leaf->left != 0 && leaf->right != 0 && node_lik_left->second.size() != node_lik_right->second.size()) {
	cerr << "Disagreement in number of characters. " << node_lik_left->second.size() << " vs. " << node_lik_right->second.size() << endl;
    }
    else {
	if (node_lik == likelihoods.end()) {
	    node_lik = likelihoods.emplace_hint(likelihoods.end(), pair<node*,vector<vector<double> > >(const_cast<node*>(leaf), vector<vector<double> >() ) );
	}
	if (leaf->left != 0) {
	    node_lik->second.resize(node_lik_left->second.size());
	    #ifdef DEBUG
	    cerr << "N chars left: " << node_lik->second.size() << endl;
	    #endif //DEBUG
	}
	else if (leaf->right != 0) {
	    node_lik->second.resize(node_lik_right->second.size());
	    #ifdef DEBUG
	    cerr << "N chars right: " << node_lik->second.size() << endl;
	    #endif //DEBUG
	}
	#ifdef DEBUG
	cerr << "N chars: " << node_lik->second.size() << endl;
	#endif //DEBUG
    }
    if (node_lik != likelihoods.end()) {
	for (unsigned int pos = 0; pos < node_lik->second.size(); ++pos) {
	    rate_model* rate_mod = partitions.get_rate_model(pos);
	    sub_model* sub_mod = partitions.get_sub_model(pos,clade_no);
	    double rate = rate_mod->get_rate_clade(clade_no);
	    node_lik->second[pos].resize(sub_mod->get_n_states());
	    if (leaf->left != 0) {
		if (node_lik_left != likelihoods.end()) {
		    branch_likelihood( node_lik->second.at(pos), node_lik_left->second.at(pos), rate_mod->get_time_mod_branch_length(leaf->left->branchlength,dist_from_root)*rate, false, *sub_mod );
		}
	    }
	    if (leaf->right != 0) {
		if (node_lik_right != likelihoods.end()) {
		    if (leaf->left != 0) branch_likelihood( node_lik->second.at(pos), node_lik_right->second.at(pos), rate_mod->get_time_mod_branch_length(leaf->right->branchlength,dist_from_root)*rate, true, *sub_mod );
		    else branch_likelihood( node_lik->second.at(pos), node_lik_right->second.at(pos), rate_mod->get_time_mod_branch_length(leaf->right->branchlength,dist_from_root)*rate, false, *sub_mod );
		}
            }
        }
    }
} 

void simpleML::branch_likelihood (  vector<double>& LHbins, vector<double>& startLHs, double branch_length, bool multiply, sub_model& model ) {
    model.set_P_matrix(branch_length);
    unsigned int n_states = model.get_n_states();
    for (unsigned int i = 0; i < n_states; ++i) {
            double pi(0.0);
            for (unsigned int j = 0; j < n_states; ++j) {
                pi += startLHs[j] * model.get_P_value(i,j);
	    }
	    if (multiply) LHbins[i] *=pi;
	    else LHbins[i] = pi;
    }
} 

void simpleML::draw_normalized_likelihood_on_nodes( node* leaf ) {
    if ( leaf!=0 ) {
	if (leaf->left != 0 || leaf->right !=0) {
	    if (leaf->left != 0) draw_normalized_likelihood_on_nodes(leaf->left);
	    if (leaf->right !=0) draw_normalized_likelihood_on_nodes(leaf->right);
	    double sum=0;
	    unsigned int char_no(0);
	    stringstream prob_label;
    	    prob_label << "[&";// P={";
	    for (vector<vector<double> >::iterator pos= likelihoods[leaf].begin(); pos != likelihoods[leaf].end(); ++pos) {
		if (pos != likelihoods[leaf].begin()) { prob_label << ","; }
		prob_label << "P={";
		for (vector<double>::iterator i= pos->begin(); i != pos->end(); ++i) sum += *i;
		for (vector<double>::iterator i= pos->begin(); i != pos->end(); ++i) {
		    if (i!=pos->begin()) prob_label << ',';
		    prob_label << *i/sum;
		    #ifdef DEBUG
		    cerr << *i << " (" << sum << ") ";
		    #endif //DEBUG
		}
		prob_label << "}";
	    }
	    prob_label << "]";
	    string label;
	    label = prob_label.str();
	    if (leaf->nodelabel!=0) label = *leaf->nodelabel+label;
	    #ifdef DEBUG
	    cerr << endl << label << endl;
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
	if (leaf->nodelabel && !leaf->nodelabel->empty()) cerr << *leaf->nodelabel << endl;
	cerr << "Simulating trait along branch. Start: " << ancestor << " end: " << trait << endl;
	#endif //DEBUG
	simdata[leaf] = trait;
	if (leaf->left != 0) simulate_chars_subtree(leaf->left, simdata, trait, model);
	if (leaf->right != 0) simulate_chars_subtree(leaf->right, simdata, trait, model);
    }
}
