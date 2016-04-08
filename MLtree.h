#include "tree.h"
#include "marth/marth.h"

using namespace std;

class MLtree : public tree {
    public:
	struct model {
	    float rate;
	    marth::square_matrix Q_matrix(n_states);
	}
	struct MLnode {
	    float rate;
	    double likelihoods[n_states]
	};
	struct char_set {
	    int start;
	    int end;
	}
	void read_states();
	int add_model(double* matrix, float rate);
	void calculate_likelihood();
    private:
	// variables
	const int n_patterns;
	const int n_states;
	vector<model> models;
	int* char_per_pattern;
	// functions
	void calc_prob_nodes(node* leaf, model* MLmodel);
}
