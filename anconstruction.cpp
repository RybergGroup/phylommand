/********************************************************************
Copyright (C) 2011 Martin Ryberg

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
#include <map>
#include <set>
#include <nlopt.hpp>
#include "simpleML.h"
#include <cfloat>
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
    char method = 'o';
    char print_tree = 'x';
    //char tree_source = 's';
    unsigned int n_states = 0;
    unsigned int n_parameters = 0;
    //int n_states
    vector<double> model_parameters;
    set<int> fixed_parameters;
    double rate_mod = 1.0;
    double cut_off = 0.0;
    vector<unsigned int> model_specifications;
    string file_name;
    ifstream input_file;
    istream* input_stream = &std::cin;
    ifstream char_input_file;
    string char_file_name;
    bool optimize_param = true;
    // Parsing arguments, present short options: -m, -o, -l, -d, -g, -i, -w, -x, -z, -h, -c, -n, -t, -u, -b, -r, -s, -p
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
	    if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--model")) {
		while (i < argc-1 && argv[i+1][0] != '-') {
		    model_specifications.push_back(atoi(argv[++i]));
		}
	    }
	    else if (!strcmp(argv[i],"-x") || !strcmp(argv[i],"--fixed")) {	
		while (i < argc-1 && argv[i+1][0] != '-') {
		    fixed_parameters.insert(atoi(argv[++i]));
		}
	    }
	    else if (!strcmp(argv[i],"-p") || !strcmp(argv[i],"--parameters")) {
                while (i < argc-1 && argv[i+1][0] != '-') {
                    model_parameters.push_back(atof(argv[++i]));
                }
            }
	    else if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--char_file")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    char_file_name = argv[++i];
                else {
                    std::cerr << "-c/--char_file require a file name as next argument." << endl;
                    return 1;
                }
            }
	    else if (!strcmp(argv[i],"-n") || !strcmp(argv[i],"--no_optim")) optimize_param = false;
	    else if (!strcmp(argv[i],"-e") || !strcmp(argv[i],"--method")) {
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
                    std::cerr << "-e/--method require the name of the method as next argument." << endl;
                    return 1;
                }
	    }
	    else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--rate_mod")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    rate_mod = atof(argv[++i]);
                else {
                    std::cerr << "-r/--rate_mod require a number as next argument." << endl;
                    return 1;
                }
            }
	    else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--time")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    cut_off = atof(argv[++i]);
                else {
                    std::cerr << "-t/--time require a number as next argument." << endl;
                    return 1;
                }
            }
            else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
                help();
                return 0;
            }
	    else if (!strcmp(argv[i],"-f") || !strcmp(argv[i],"--file")) {
		if ( i < argc-1 && argv[i+1][0] != '-' )
		    file_name = argv[++i];
		else {
		    std::cerr << "-f/--file require a file name as next argument" << endl;
		    return 1;
		}
	    }
	    else if (i == argc-1 && argv[i][0] != '-' ) {
		file_name = argv[i];
	    }
            else {
                std::cerr << "Unrecognized argument " << argv[i] << ". Quitting quietly." << endl;
                return 1;
            }
        }
    }
    if (!file_name.empty()) {
	input_file.open(file_name.c_str(),std::ifstream::in);
	if (input_file.good())
	    input_stream = &input_file;
	else {
	    cerr << "Could not open file: " << file_name << endl;
	    return 1;
	}
    }
    // set characters
    map<string,int> characters;
    if (!char_file_name.empty()) {
        char_input_file.open(char_file_name.c_str(),std::ifstream::in);
        if (char_input_file.good()) {
	    string word;
	    string taxon;
	    string character;
	    set<unsigned int> character_states;
            while (char_input_file >> word) {
		if (taxon.empty()) taxon = word;
		else character = word;
		if (!character.empty()) {
		    unsigned int state = atoi(character.c_str());
		    characters[taxon] = state;
		    character_states.insert(state);
		    taxon.clear();
		    character.clear();
		}
	    }
	    n_states = character_states.size();
	    for (unsigned int i=0; i < n_states; ++i)
		if (character_states.find(i) == character_states.end()) cerr << "No taxa with character " << i << "!" << endl;
	}
        else {
            cerr << "Could not open file: " << char_file_name << endl;
            return 1;
        }
	if (char_input_file.is_open()) char_input_file.close();
    }
    else {
	cerr << "Could not find file name." << endl;
	return 1;
    }
    // Set number of parameters in rate matrix
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
    else if (n_parameters<0) {
	std::cerr << "Error, negative number of free parameters!" << endl;
	return 1;
    }
    #if DEBUG
	std::cerr << "2. Number of fixed parameters: " << fixed_parameters.size() << endl;
	std::cerr << "2. Fixed parameters: ";
	for (set<int>::iterator i=fixed_parameters.begin(); i!=fixed_parameters.end(); ++i) std::cerr << *i << ' ';
	std::cerr << endl;
	std::cerr << "2. Number of parameters: " << n_parameters << endl;
    #endif //DEBUG
    while (1) {
	simpleML tree;
	// Get tree
	tree.tree_file_parser( *input_stream );
	if (tree.n_tips() < 2) break;
        tree.init(n_states);
	for (map<string,int>::iterator i= characters.begin(); i!=characters.end();++i) {
	    double state_prob[n_states];
	    for (unsigned int j=0; j<n_states; ++j) state_prob[j] = 0.0;
	    state_prob[i->second]=1.0;
	    tree.set_char(i->first,state_prob);
	}
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
	    if (print_tree == 'x') std::cout << '[';
	    if (method == 'o') {
		maximize.set_max_objective(opt_function,&data);
		cout << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
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
		cout << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
		//double LogL;
		nlopt::result result = maximize.optimize(variable_values,LogL);
		if (result < 0) {
		    cerr << "Failure when optimizing!!!" << endl;
		    return 1;
		}
		cout << "Time from root to rate change: " << model_parameters[model_parameters.size()-2] << endl;
		cout << "Rate multiplier: " << model_parameters[model_parameters.size()-1] << endl;
	    }
	    else {
		cerr << "Unrecognized method. Nothing to do." << endl;
		return 1;
	    }
	    cout << "Parameter values:" << endl;
	    #if DEBUG
    	    for (unsigned int i= 0; i < model_parameters.size(); ++i) {
		cout << i << ": " << model_parameters[i] << ' ';
	    }
	    cout << endl;
	    #endif /*DEBUG*/
	    cout << "Q matrix:" << endl;
	    tree.print_Q_matrix();
	    cout << endl;
	    cout << "Log likelihood: " << LogL << endl;
	}
	else {
	    if (method == 'o') {
		tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
		cout << "Log likelihood= " << tree.calculate_log_likelihood() << endl;
	    }
	    else if (method == 't') {
                tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
                cout << "Log likelihood= " << tree.calculate_likelihood_rate_change_in_time(cut_off,rate_mod) << endl;
		cout << "Time from root: " << cut_off << " Rate modifier: " << rate_mod << endl;
            }
	    else {
                cerr << "Unrecognized method. Nothing to do." << endl;
                return 1;
            }
	    tree.print_Q_matrix();
    	    cout << endl;
	}
	if (print_tree == 'w') {
	    cout << "Tree:" << endl;
	    tree.draw_normalized_likelihood_on_nodes();
	    tree.print_newick();
	}
	else if (print_tree == 'x') {
	    std::cout << ']' << std::endl;
	    tree.draw_normalized_likelihood_on_nodes();
	    tree.print_nexus();
	}
    }
    // Close file stream if open
    if (input_file.is_open()) input_file.close();
}

void help () {
    std::cout << "anconstruction is a command line program for calculating likelihood on a trees." << endl;
    std::cout << "The program take a tree in newick format as indata through standard in." << endl;
    std::cout << "(c) Martin Ryberg 2015." << endl << endl;
    std::cout << "Usage:" << endl << "anconstruction [arguments] < file.tree" << endl;
    std::cout << "anconstruction [arguments] file.tree" << endl << endl;
    std::cout << "For the second alternative you need to be careful so treebender does not interpret the filename as an extra argument to a switch. If this happen treebender will expect" <<endl;
    std::cout << "input from standard in and it will appear as nothing is happening. This can be avoided by giving the filename after the switch --file/-f (see below)." <<endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--char_file/-c [file name]                                 give file name of character file, e.g. -c character.file." << endl;
    std::cout << "--file/-f [file name]                                      give name of tree file, e.g. -f file.tree." << endl;
    std::cout << "--help / -h                                                print this help." << endl;
    std::cout << "--method/-e [name of method]                               set which method to use. The method name should be next argument, e.g. -e const. Default: const." << endl;
    std::cout << "    const                                                      maximize likelihood for parameters with constant rate through the tree." << endl;
    std::cout << "    time                                                       maximize likelihood for given parameters with rate shift at one point in time." << endl;
    std::cout << "--model/-m [space separated string of integer numbers]     give the model by numbering the rate parameters for different transition, e.g. -m 0,1,0,2,1,2" << endl;
    std::cout << "--no_optim/-n                                              calculate likelihood for given parameters. No optimization." << endl;
    std::cout << "--parameters/-p [space separated string of real numbers]   give corresponding parameter values for parameters. If optimizing these will be starting values, e.g. -p 0.1 0.01 0.05" << endl;
    std::cout << "--rate_mod/-r [real number]                                give modifier for rate compared to rate at root, e.g. -r 0.5. Default: 1.0." << endl;
    std::cout << "--time/-t [real number]                                    give branch length distance from root where change in rate occur, e.g. -t 10. Default: 0" << endl;
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
