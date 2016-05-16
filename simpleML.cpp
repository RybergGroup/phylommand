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
	//leaf->other = 0;
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

void simpleML::set_Q_matrix ( const double* values ) {
    unsigned int parameters[(n_states*n_states)-n_states];
    for (unsigned int i=0; i<(n_states*n_states)-n_states; ++i) parameters[i]=i;
    set_Q_matrix (&parameters[0], values);
}

void simpleML::set_Q_matrix ( const unsigned int* parameters, const double* values ) {
    for (unsigned char i=0; i < n_states; ++i) {
	for (unsigned char j=0; j < n_states; ++j) { 
	    if (i == j) {
		double value = 0;
		for (unsigned int k = (i*n_states)-i; k < (i+1)*(n_states-1); ++k) value -= values[parameters[k]];
		Q_matrix.set_value(i,j,value);
	    }
	    else {
		if (j>i) Q_matrix.set_value(i,j,values[parameters[i*n_states+j-i-1]]);
		else Q_matrix.set_value(i,j,values[parameters[i*n_states+j-i]]);
	    }
	}
    }
}

void simpleML::calculate_likelihood (const node* leaf) {
    if (leaf->left != 0) {
	calculate_likelihood (leaf->left);
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),leaf->left->branchlength,false);
    }
    if (leaf->right != 0) {
	calculate_likelihood (leaf->right);
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength,true);
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),leaf->right->branchlength,false);
    }
}

void simpleML::calculate_likelihood_rate_change_in_time(const node* leaf, const double dist_from_root, const double cut_off, double rate) {
    if (leaf->left != 0) {
	calculate_likelihood_rate_change_in_time(leaf->left, dist_from_root+leaf->left->branchlength,cut_off,rate);
	double length;
	if (dist_from_root+leaf->left->branchlength > cut_off) {
	    if (dist_from_root>cut_off)
		length = leaf->left->branchlength * rate;
	    else 
		length = (leaf->left->branchlength - (cut_off-dist_from_root))*rate + cut_off-dist_from_root;
	}
	else
	    length = leaf->left->branchlength;
	branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->left)),length,false);
    }
    if (leaf->right != 0) {
	calculate_likelihood_rate_change_in_time(leaf->right, dist_from_root+leaf->right->branchlength,cut_off,rate);
	double length;
	if (dist_from_root+leaf->right->branchlength > cut_off) {
            if (dist_from_root>cut_off)
                length = leaf->right->branchlength * rate;
            else 
                length = (leaf->right->branchlength - (cut_off-dist_from_root))*rate + cut_off-dist_from_root;
        }
	else
    	    length = leaf->right->branchlength;
	if (leaf->left != 0) branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),length,true);
	else branch_likelihood(likelihoods.find(const_cast<node*>(leaf)),likelihoods.find(const_cast<node*>(leaf->right)),length,false);
    }
}

void simpleML::branch_likelihood ( map<node*,vector<double> >::iterator LHbins, map<node*,vector<double> >::iterator startLHs, double branch_length, bool multiply) {
    marth::square_matrix P_matrix;
    Q_matrix.exponential(&P_matrix, branch_length, 20);
    for (unsigned int i = 0; i < n_states; ++i) {
            double pi = 0;
            for (unsigned int j = 0; j < n_states; ++j)
                pi += startLHs->second[j] * P_matrix.get_value(i,j);
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
