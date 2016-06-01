/********************************************************************
Copyright (C) 2013 Martin Ryberg

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

contact: kryberg@utk.edu
*********************************************************************/

#include <iostream>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <cfloat>
#include <nlopt.hpp>
#include "tree.h"
#include "nj_tree.h"
#include "character_vector.h"
#include "matrix_parser.h"
#include "simpleML.h"

using namespace std;

void help ();

double opt_function(const std::vector<double> &x, std::vector<double> &grad, void* data);
double opt_rate_in_time(const std::vector<double> &x, std::vector<double> &grad, void* data);
void change_non_fixed(const std::vector<double> &x, std::vector<double>* values, const std::set<int>* fixed);

struct tree_modelspec_struct {
    simpleML* tree;
    const vector<unsigned int>* specifications;
    vector<double>* values;
    const set<int>* fixed;
};

int main (int argc, char *argv []) {
    #if DEBUG
    cerr << "Debugging version!!!" << endl;
    #endif //DEBUG
    bool lables = true;
    char method = 'p';
    bool random = false;
    bool print_br_length(true);
    bool print_state_on_nodelable(false);
    string tree_file_name;
    ifstream tree_file;
    istream* tree_stream = &std::cin;
    string data_file_name;
    ifstream data_file;
    istream* data_stream = &std::cin;
    string alphabet_file_name;
    ///// Variables for ancon
    bool optimize_param = true;
    vector<unsigned int> model_specifications;
    double cut_off = 0.0;
    double rate_mod = 1.0;
    set<int> fixed_parameters;
    vector<double> model_parameters;
    char print_tree = 'x';
    //////////////////////////
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-L") || !strcmp(argv[i],"--no_lable")) lables = false;
            else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--data_file")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    data_file_name = argv[++i];
                else {
                    std::cerr << "-d/--data_file require a file name as next argument" << endl;
                    return 1;
                }
	    }
            else if (!strcmp(argv[i],"-n") || !strcmp(argv[i],"--neighbour_joining")) {
		method = 'n';
	    }
            else if (!strcmp(argv[i],"-p") || !strcmp(argv[i],"--parsimony")) {
		method = 'p';
	    }
            //else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--likelihood")) {
	//	method = 'l';
	  //  }
	    else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--step_wise")) {
		method = 's';
	    }
            else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
                help();
                return 0;
            }
	    else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--tree_file")) {
		//if (indata == '0') indata = 'd';
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    tree_file_name = argv[++i];
                else {
                    std::cerr << "-t/--tree_file require a file name as next argument" << endl;
                    return 1;
                }
            }
	    else if (!strcmp(argv[i],"-a") || !strcmp(argv[i],"--alphabet_file")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    alphabet_file_name = argv[++i];
                else {
                    std::cerr << "-a/--alphabet_file require a file name as next argument" << endl;
                    return 1;
                }
            }
	    else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--random")) {
                random = true;
            }
            else if (!strcmp(argv[i],"-0") || !strcmp(argv[i],"--no_branch_length")) print_br_length = false;
///////////// From ancon
            else if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--model")) {
                while (i < argc-1 && argv[i+1][0] != '-') {
                    model_specifications.push_back(atoi(argv[++i]));
                }
            }
            else if (!strcmp(argv[i],"-x") || !strcmp(argv[i],"--fixed")) {
                while (i < argc-1 && argv[i+1][0] != '-') {
                    fixed_parameters.insert(atoi(argv[++i]));
                }
            }
            else if (!strcmp(argv[i],"-P") || !strcmp(argv[i],"--parameters")) {
                while (i < argc-1 && argv[i+1][0] != '-') {
                    model_parameters.push_back(atof(argv[++i]));
                }
            }
            else if (!strcmp(argv[i],"-N") || !strcmp(argv[i],"--no_optim")) optimize_param = false;
            else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--likelihood")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) {
                    string method_name = argv[++i];
                    if (!method_name.compare("const")) method = 'o';
                    else if (!method_name.compare("time")) method='t';
                    else {
                        std::cerr << "Do not recognize method \"" << method_name << "\"." << endl;
                        return 1;
                    }
                }
                else {
		    method = 'o';
                    //std::cerr << "-e/--method require the name of the method as next argument." << endl;
                    //return 1;
                }
            }
            else if (!strcmp(argv[i],"-R") || !strcmp(argv[i],"--rate_mod")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    rate_mod = atof(argv[++i]);
                else {
                    std::cerr << "-R/--rate_mod require a number as next argument." << endl;
                    return 1;
                }
            }
            else if (!strcmp(argv[i],"-T") || !strcmp(argv[i],"--time")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    cut_off = atof(argv[++i]);
                else {
                    std::cerr << "-T/--time require a number as next argument." << endl;
                    return 1;
                }
            }
	    else if (!strcmp(argv[i],"--get_state_at_nodes")) {
		print_state_on_nodelable = true;
	    }
////////////////////////
            else if (i == argc-1 && argv[i][0] != '-' ) {
		if (data_file_name.empty()) data_file_name = argv[i];
                else if (tree_file_name.empty()) tree_file_name = argv[i];
		else {
		    std::cerr << "Have no use for argument " << argv[i] << ". Already have file names for tree and data." << std::endl;
		    return 1;
		}
            }
            else {
                std::cerr << "Unrecognized argument " << argv[i] << ". Quitting quietly." << endl;
                return 1;
            }
        }
    }
    if (!tree_file_name.empty()) {
        tree_file.open(tree_file_name.c_str(),std::ifstream::in);
        if (tree_file.good())
            tree_stream = &tree_file;
        else {
            cerr << "Could not open file: " << tree_file_name << endl;
            return 1;
        }
    }
    if (!data_file_name.empty()) {
        data_file.open(data_file_name.c_str(),std::ifstream::in);
        if (data_file.good())
            data_stream = &data_file;
        else {
            cerr << "Could not open file: " << data_file_name << endl;
            return 1;
        }
    }
    if (method == 'n') {
	njtree tree;
	tree.read_distance_matrix(*data_stream, lables);
	//tree.print_node_and_distance_array();
	tree.build_nj_tree();
	tree.print_newick(print_br_length);
	return 0;
    }
    if (method == 'p' || method == 's' || method == 'o' || method == 't') {
	vector<character_vector> characters;
////////////////////
	map<char, bitset<SIZE> > alphabet;
	alphabet::set_alphabet_binary(alphabet);
	if (!alphabet_file_name.empty()) {
	    ifstream alphabet_file;
	    alphabet_file.open(alphabet_file_name.c_str(),std::ifstream::in);
	    alphabet.clear();
	    alphabet_parser parser(alphabet_file,alphabet);
	    parser.pars();
	}
	#ifdef DEBUG
	cerr << "Number of different characters: " << alphabet.size() << endl;
	for (map<char,bitset<SIZE> >::const_iterator i = alphabet.begin(); i != alphabet.end(); ++i)
	    cerr << i->first << " - " << i->second << endl;
	#endif //DEBUG
	partitions regions;
	regions.add_alphabet("first",alphabet);
	regions.add_partition(0,0,"default","first");
	matrix_parser data_parser(*data_stream, characters, regions);
	data_parser.pars();
	#ifdef DEBUG
	std::cerr << "Number of taxa: " << characters.size() << endl;
	std::cerr << "Max number of characters: " << characters.begin()->max_n_char() << std::endl;
	#endif //DEBUG
        if (data_file.is_open()) data_file.close();
/////////////////////////////////
	if (method == 'p') {
	    while (1) {
		tree tree;
		tree.tree_file_parser( *tree_stream );
		if (tree.empty()) break;
		std::cout << tree.fitch_parsimony( characters, print_state_on_nodelable, print_br_length, alphabet ) << std::endl;
		if (print_br_length) tree.print_newick(print_br_length);
	    }
	}
	else if (method == 's') {
	    tree tree;
	    if (random) random_shuffle(characters.begin(),characters.end());
	    tree.stepwise_addition(characters);
	    if (print_br_length) tree.fitch_parsimony( characters, print_br_length );
	    tree.print_newick(print_br_length);
	}
	if (method == 'o' || method == 't') {
	    unsigned int n_states = 0;
	    unsigned int n_parameters = 0;
	    for (vector<character_vector>::iterator i=characters.begin(); i!=characters.end(); ++i)
		if (i->highest_char_state()+1 > n_states) n_states = i->highest_char_state()+1;
	    if (!model_specifications.empty()) {
		set<unsigned int> parameters;
		for (vector<unsigned int>::const_iterator i=model_specifications.begin(); i!=model_specifications.end(); ++i) {
		    parameters.insert(*i);
		}
		n_parameters=parameters.size();
	    }
	    else {
		for (unsigned int i=0; i < (n_states*n_states)-n_states; ++i) model_specifications.push_back(i);
		n_parameters=model_specifications.size();
	    }
// Set parameter values
	    if(!model_parameters.empty()) {
		if (model_parameters.size() < n_parameters) {
		    cerr << "The number of parameter values (" << model_parameters.size() << ") must be as many as parameters (" << n_parameters << ")." << endl;
		    return 1;
		}
		for (vector<double>::const_iterator value=model_parameters.begin(); value!=model_parameters.end(); ++value) {
		    if (*value < 0.0) {
			cerr << "All parameter values must be positive." << endl;
			return 1;
		    }
		}
	    }
	    else {
		for (unsigned int i=0; i < n_parameters; ++i) model_parameters.push_back(0.1);
	    }
	    // Adjust number of parameters according to rate change model and add extra parameters to vector
	    if (method == 't') {
		n_parameters+=2;
		model_parameters.push_back(cut_off);
		model_parameters.push_back(rate_mod);
	    }
	    #if DEBUG
		std::cerr << "1. Number of fixed parameters: " << fixed_parameters.size() << endl;
		std::cerr << "1. Number of parameters: " << n_parameters << endl;
	    #endif //DEBUG
	    n_parameters-=fixed_parameters.size();
	    if (n_parameters==0)
		optimize_param=false;


///////////////
	    while (1) {
		stringstream model_out;
		simpleML tree;
		tree.tree_file_parser( *tree_stream );
		if (tree.empty()) break;
		tree.init(n_states); // need to fix this
		for (vector<character_vector>::iterator i=characters.begin(); i!=characters.end(); ++i)
		    tree.set_char(i->get_taxon(),i->get_character(0));
////////////////////////////////////////////
		if (optimize_param) {
		    nlopt::opt maximize(nlopt::LN_NELDERMEAD, n_parameters);
		    tree_modelspec_struct data = {&tree, &model_specifications, &model_parameters, &fixed_parameters};
		    vector<double> variable_values;
		    for (int i=0; i<model_parameters.size(); ++i) {
			if (fixed_parameters.find(i)==fixed_parameters.end()) variable_values.push_back(model_parameters[i]);
		    }
		    vector<double> lower_bounds;
		    for (unsigned int i=0; i < n_parameters; ++i) lower_bounds.push_back(0.0);
		    maximize.set_lower_bounds(lower_bounds);
		    double LogL;
		    if (print_tree == 'x') model_out << '[';
		    if (method == 'o') {
			maximize.set_max_objective(opt_function,&data);
			model_out << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			nlopt::result result = maximize.optimize(variable_values,LogL);
			if (result < 0) {
			    cerr << "Failure when optimizing!!!" << endl;
			    return 1;
			}
		    }
		    else if (method == 't') {
			vector<double> upper_bounds;
			for (unsigned int i=0; i < n_parameters; ++i) {
			    if (i == n_parameters-2)
				upper_bounds.push_back(tree.longest_to_tip());
			    else
				upper_bounds.push_back(DBL_MAX);
			}
			maximize.set_upper_bounds(upper_bounds);
			maximize.set_max_objective(opt_rate_in_time,&data);
			model_out << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			//double LogL;
			nlopt::result result = maximize.optimize(variable_values,LogL);
			if (result < 0) {
			    cerr << "Failure when optimizing!!!" << endl;
			    return 1;
			}
			model_out << "Time from root to rate change: " << model_parameters[model_parameters.size()-2] << endl;
			model_out << "Rate multiplier: " << model_parameters[model_parameters.size()-1] << endl;
		    }
		    else {
			cerr << "Unrecognized method. Nothing to do." << endl;
			return 1;
		    }
		    model_out << "Parameter values:" << endl;
		    #if DEBUG
		    for (unsigned int i= 0; i < model_parameters.size(); ++i) {
			model_out << i << ": " << model_parameters[i] << ' ';
		    }
		    model_out << endl;
		    #endif /*DEBUG*/
		    model_out << "Q matrix:" << endl;
		    tree.print_Q_matrix( model_out );
		    model_out << endl;
		    model_out << "Log likelihood: " << LogL << endl;
/////////////////////
//		    cerr << "Not implemented yet" << endl;
		}
		/////////////////////////////
		else {
		    if (method == 'o') {
			tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
			model_out << "Log likelihood= " << tree.calculate_log_likelihood() << endl;
		    }
		    else if (method == 't') {
			tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
			model_out << "Log likelihood= " << tree.calculate_likelihood_rate_change_in_time(cut_off,rate_mod) << endl;
			model_out << "Time from root: " << cut_off << " Rate modifier: " << rate_mod << endl;
		    }
		    else {
			cerr << "Unrecognized method. Nothing to do." << endl;
			return 1;
		    }
		    tree.print_Q_matrix( model_out );
		    model_out << endl;
		}
		if (print_tree == 'w') {
		    // cout << "Tree:" << endl;
		    if (print_state_on_nodelable) tree.draw_normalized_likelihood_on_nodes();
		    tree.print_newick();
		}
		else if (print_tree == 'x') {
		    model_out << ']' << std::endl;
		    if (print_state_on_nodelable) tree.draw_normalized_likelihood_on_nodes();
		    tree.print_nexus();
		}
		cout << model_out.str();
/////////////////////////////
	    }
	    //std::cout << "Sorry, likelihood is not available yet." << std::endl;
	    //return 0;
	}
    }
}

void help () {
    std::cout << "Treeator is a command line program to construct trees." << endl;
    std::cout << "The program take either a left triangular similarity" << endl;
    std::cout << "matrix (neighbour joining) or a data matrix of relaxed" << endl;
    std::cout << "phylip format (not interleaved; parsimony/maximum likelihood)" << endl;
    std::cout << " as input through standard in/ last argument/ as given below." << endl;
    std::cout << "(c) Martin Ryberg 2015." << endl << endl;
    std::cout << "Usage:" << endl << "treeator [arguments] < data_file.txt" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--alphabet_file / -a                                       give file with character alphabet." << std::endl;
    std::cout << "--data_file / -d [file name]                               give the data file."<< std::endl;
    std::cout << "--fixed/-x [parameter name]                                give name of parameter to fix. This argument can be repeated." << endl;
    std::cout << "--get_state_at_nodes                                       will give the states at nodes as comments (readable in FigTree)." << endl;
    std::cout << "--help / -h                                                print this help." << endl;
    std::cout << "--likelihood / -l                                          calculate likelihood for data given tree. Either with constant rate throug time (const) or" << endl;
    std::cout << "                                                               with rate changing (multiplied by a variable) at a certain time point (time)." << endl;
    std::cout << "--model/-m [space separated string of integer numbers]     give the model by numbering the rate parameters for different transition, e.g. -m 0,1,0,2,1,2" << endl;
    std::cout << "--neighbour_joining / -n                                   compute neighbour joining tree for given data. The data should be" << std::endl;
    std::cout << "                                                           a left triangular similarity matrix." << std::endl;
    std::cout << "--no_branch_length / -0                                    do not print branch lengths and do not calculate branch lengths for" << std::endl;
    std::cout << "                                                               parsimony trees" << endl;
    std::cout << "--no_lable / -l                                            will tell treeator that there are no taxon labels in the matrix." << endl;
    std::cout << "--no_optim/-n                                              calculate likelihood for given parameters. No optimization." << endl;
    std::cout << "--parameters/-P [space separated string of real numbers]   give corresponding parameter values for parameters. If optimizing these will be starting values," << endl;
    std::cout << "                                                               e.g. -p 0.1 0.01 0.05" << endl;
    std::cout << "--parsimony / -p                                           calculate parsimony score for given tree and data." << std::endl;
    std::cout << "--rate_mod/-R [real number]                                give modifier for rate compared to rate at root, e.g. -r 0.5. Default: 1.0." << endl;
    std::cout << "--random / -r                                              do stepwise addition in random order." << endl;
    std::cout << "--time/-T [real number]                                    give branch length distance from root where change in rate occur, e.g. -t 10. Default: 0" << endl;
    std::cout << "--tree_file / -t [file name]                               give tree file name" << std::endl;
    std::cout << "--step_wise / -s                                           do parsimony stepwise addition." << std::endl;
    std::cout << endl;
}

double opt_function(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    change_non_fixed(x,static_cast<tree_modelspec_struct*>(data)->values,static_cast<tree_modelspec_struct*>(data)->fixed);
    double LogLH = static_cast<tree_modelspec_struct*>(data)->tree->calculate_log_likelihood(&(static_cast<tree_modelspec_struct*>(data)->specifications->front()),&(static_cast<tree_modelspec_struct*>(data)->values->at(0)));
    #if DEBUG
        std::cerr << "Likelihood: " << LogLH << endl;
    #endif //DEBUG
    return LogLH;
}

double opt_rate_in_time(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    change_non_fixed(x,static_cast<tree_modelspec_struct*>(data)->values,static_cast<tree_modelspec_struct*>(data)->fixed);
    if (static_cast<tree_modelspec_struct*>(data)->values->size()<3) return 0.0;
    double LogLH =  static_cast<tree_modelspec_struct*>(data)->tree->calculate_likelihood_rate_change_in_time(&(static_cast<tree_modelspec_struct*>(data)->specifications->front()),&(static_cast<tree_modelspec_struct*>(data)->values->at(0)),static_cast<tree_modelspec_struct*>(data)->values->at(static_cast<tree_modelspec_struct*>(data)->values->size()-2),static_cast<tree_modelspec_struct*>(data)->values->at(static_cast<tree_modelspec_struct*>(data)->values->size()-1));
    #if DEBUG
        std::cerr << "Likelihood: " << LogLH << endl;
    #endif //DEBUG
    return LogLH;
}

void change_non_fixed(const std::vector<double> &x, std::vector<double>* values, const std::set<int>* fixed) {
    #if DEBUG
        std::cerr << "No. fixed "<< fixed->size();
        std::cerr << " (";
        for (set<int>::iterator i=fixed->begin(); i!=fixed->end(); ++i) std::cerr << *i << ' ';
        std::cerr << ')' << endl;
        std::cerr << "No. parameters " << values->size() << endl;
        std::cerr << "Values to assign: ";
        for (vector<double>::const_iterator i=x.begin(); i!=x.end(); ++i) std::cerr << *i << ' ';
        std::cerr << endl;
    #endif //DEBUG
    int z=0;
    int j=0;
    for (std::vector<double>::iterator i=values->begin(); i != values->end(); ++i) {
        #if DEBUG
            std::cout << "Param " << z << " has value " << *i << endl;
        #endif //DEBUG
        if (fixed->find(z)==fixed->end()) {
            #if DEBUG
                cout << "Assigning " << x[j] << '(' << j << ") to parameter with value " << *i << '(' << z << ')'<< endl;
            #endif //DEBUG
            *i = x[j];
            ++j;
        }
        ++z;
    }
}

