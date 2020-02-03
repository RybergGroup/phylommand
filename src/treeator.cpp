/********************************************************************
Copyright (C) 2019 Martin Ryberg

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
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <map>
#include <cfloat>
#ifdef NLOPT
#include <nlopt.hpp>
#endif //NLOPT
#include "tree.h"
#include "nj_tree.h"
#include "character_vector.h"
#include "file_parser.h"
#include "argv_parser.h"
#include "matrix_parser.h"
#include "simpleML.h"
#include "sub_model.h"
#include "rate_model.h"

using namespace std;

void help ();

double opt_function(const std::vector<double> &x, std::vector<double> &grad, void* data);
double opt_rate_in_time(const std::vector<double> &x, std::vector<double> &grad, void* data);
double opt_rate_for_clades_function(const std::vector<double> &x, std::vector<double> &grad, void* data);
//void change_non_fixed(const std::vector<double> &x, std::vector<double>* values, const std::set<unsigned int>* fixed);

struct tree_modelspec_struct {
    simpleML* tree;
    sub_model* model;
    const set<unsigned int>* fixed;
    const bool fixed_freq;
    const unsigned int freq_start;
    //const vector<unsigned int>* specifications; // which parameter represents which rate i.e. 0,0,1,1,2,2 represents same rates going from 0->1 as 0->2, and 1->0 as 1->2, and so on.
    rate_model* rate_mod;
    const unsigned int extra_start;
    const set<unsigned int>* fixed_extra;
    //const vector<unsigned int>* extra_specs;
    //bool tr;
};

bool change_non_fixed(const std::vector<double> &x, tree_modelspec_struct* data);

int main (int argc, char *argv []) {
    #if DEBUG
    cerr << "Debugging version!!!" << endl;
    #endif //DEBUG
    bool labels = true;
    char method = 'p';
    bool print_br_length(true);
    bool calc_br_length(true);
    bool print_state_on_nodelabel(false);
    char print_tree = 'w';
    bool quiet(true);
    bool tr(false);
    bool non_as_uncertain(true);
    unsigned int n_char = 1;
    string tree_file_name;
    ifstream tree_file;
    file_parser tree_input(&cin);
    string data_file_name;
    ifstream data_file;
    file_parser data_input(&cin);
    string alphabet_file_name;
    string tree_file_format;
    string data_file_format;
    map<string,string> taxa_trans;
    vector<string> taxon_sets;
    ///// Variables for ancon
    #ifdef NLOPT
    bool optimize_param = true;
    #else
    bool optimize_param = false;
    #endif //NLOPT
    vector<unsigned int> model_specifications;
    //vector<unsigned int> rate_mod_specs;
    //double cut_off = 0.0;
    rate_model rate_mod;
    set<unsigned int> fixed_parameters;
    set<unsigned int> fixed_extra_parameters;
    bool fixed_freq(false);
    vector<double> model_parameters;
    vector<double> char_freqs;
    ///////////
    bool tree_same_as_data_when_nexus(false);
    /// Reading arguments ///
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-L") || !strcmp(argv[i],"--no_label")) labels = false;
            else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--data_file")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    if (data_file_name.empty())
			data_file_name = argv[i];
		    else
			cerr << "Data file already given (" << data_file_name << "). Will ignore " << argv[i] << "." << endl;
		}
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
	    else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--step_wise")) {
		method = 's';
	    }
	    else if (!strcmp(argv[i],"--simulate")) {
		method = 'C';
		if ( i < argc-1 && argv[i+1][0] != '-' )
		    n_char = atoi(argv[++i]);
	    }
    	    else if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose")) quiet = false;
            else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
                help();
                return 0;
            }
	    else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--tree_file")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    if (tree_file_name.empty())
			tree_file_name = argv[i];
		    else
			cerr << "Tree file already given (" << tree_file_name << "). Will ignore " << argv[i] << '.' << endl;
		}
                else {
                    std::cerr << "-t/--tree_file require a file name as next argument" << endl;
                    return 1;
                }
            }
	    else if (!strcmp(argv[i],"-a") || !strcmp(argv[i],"--alphabet_file")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    alphabet_file_name = argv[++i];
                else {
                    std::cerr << "-a/--alphabet_file require a file name, or dna, protein or binary for dna, amino acid respectively binary alphabet, as next argument" << endl;
                    return 1;
                }
		#ifdef DEBUG
		cerr << "Got which file holds the alphabet: " << alphabet_file_name << endl;
		#endif //DEBUG
            }
	    else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--random")) {
                method = 'r';
            }
            else if (!strcmp(argv[i],"-0") || !strcmp(argv[i],"--no_branch_length")) print_br_length = false;
	    else if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--do_not_calc_br_length")) calc_br_length = false;
///////////// From ancon
            else if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--model")) {
                if (i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    string number;
		    unsigned int j(0);
		    while (true) {
			if ((argv[i][j] == '\0' || argv[i][j] == ',') && !number.empty()) {
			    model_specifications.push_back(atoi(number.c_str()));
			    number.clear();
			}
			else number += argv[i][j];
			if (argv[i][j] == '\0') break;
			++j;
		    }
		}
		else {
		    cerr << "-m / --model require a comma separated integer string as next argument.";
		    return 1;
		}
            }
            else if (!strcmp(argv[i],"-e") || !strcmp(argv[i],"--fixed")) {
                if (i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ',', arguments);
		    for (vector<string>::const_iterator i=arguments.begin(); i != arguments.end(); ++i)
			fixed_parameters.insert(atoi(i->c_str()));
		}
            }
	    else if (!strcmp(argv[i],"--fix_freq")) {
		fixed_freq = true;
	    }
            else if (!strcmp(argv[i],"--fixed_extras")) {
                if (i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ',', arguments);
		    for (vector<string>::const_iterator i=arguments.begin(); i != arguments.end(); ++i)
			fixed_extra_parameters.insert(atoi(i->c_str()));
		}
            }
            else if (!strcmp(argv[i],"-P") || !strcmp(argv[i],"--parameters")) {
                if (i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> numbers;
		    argv_parser::pars_sub_args(argv[i], ',', numbers);
		    for (vector<string>::const_iterator i=numbers.begin(); i != numbers.end(); ++i)
	    		model_parameters.push_back(atof(i->c_str()));
		}
		else {
		    cerr << "-P/--parameters require a comma separated real number string as next argument.";
		    return 1;
		}
            }
	    else if (!strcmp(argv[i],"-F") || !strcmp(argv[i],"--frequencies")) {
		if (i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> numbers;
		    argv_parser::pars_sub_args(argv[i], ',', numbers);
		    for (vector<string>::iterator sub = numbers.begin(); sub != numbers.end(); ++sub) {
                        char_freqs.push_back(atof(sub->c_str()));
		    }
		}
                else {
                    cerr << "-F/--frequencies require a comma separated real number string as next argument.";
                    return 1;
                }
		if (method != 'A' && method != 't' && method !='S' && method != 'C') method = 'o';
	    }
	    else if (!strcmp(argv[i],"--time_reversible") || !strcmp(argv[i],"--tr")) tr = true; 
            else if (!strcmp(argv[i],"-N") || !strcmp(argv[i],"--no_optim")) optimize_param = false;
            else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--likelihood")) {
		    if (method != 't' && method != 'A' && method !='S') method = 'o';
            }
	    else if (!strcmp(argv[i],"-U") || !strcmp(argv[i],"--rate_change_at_time")) {
                method = 't';
                if ( i < argc-1 && argv[i+1][0] != '-') {
                    ++i;
                    vector<string> arguments;
                    argv_parser::pars_sub_args(argv[i], ':', arguments );
                    rate_mod.set_time_rate(0, atof(arguments[0].c_str()));
                    if (!rate_mod.add_rate_for_time(0,atof(arguments[1].c_str()))) {
                        cerr << "Could not add rate " << atof(arguments[0].c_str()) << " for cut off " << atof(arguments[1].c_str()) << endl;
                        return 1;
                    }
                }
            }
	    else if (!strcmp(argv[i],"-S") || !strcmp(argv[i],"--skyline_rates")) {
		method = 'S';
		if ( i < argc-1 && argv[i+1][0] != '-') {
                    ++i;
                    if (!rate_mod.pars_rates_and_times_arg(argv[i])) {
                        cerr << "Unable to pars \"" << argv[i] << "\" for rates" << endl;
                        return 1;
                    }
		}
		else { cerr << "-S/--skyline_rates require a list of time cut offs and rate multipliers as second argument." << endl; return 1; }	
	    }
	    else if (!strcmp(argv[i],"-V") || !strcmp(argv[i],"-A") || !strcmp(argv[i],"--clade_rates")) {
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
                    if (!rate_mod.pars_clade_rates_arg(argv[i],taxon_sets)) {
                        cerr << "Unable to pars \"" << argv[i] << "\" for rates" << endl;
                        return 1;
                    }
		}
		else {
		    cerr << "-V / --clade_rates require at least one set of taxa given as a comma separated string as next argument." << endl;
		    return 1;
		}
		method = 'A';
	    }
	    else if (!strcmp(argv[i],"--get_state_at_nodes")) {
		print_state_on_nodelabel = true;
	    }
            else if (!strcmp(argv[i],"--format")) {
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    string first; string second;
		    string* taking_it = &first;
		    for (unsigned int j=0; argv[i][j] != '\0'; ++j) {
			if (argv[i][j] == ',') taking_it = &second;
			else *taking_it += argv[i][j];
		    }
		    if (!first.compare("nexus")) data_file_format = "nexus";
		    else if (!first.compare("fasta")) data_file_format = "fasta";
		    else if (!first.compare("phylip")) data_file_format = "phylip";
		    else {
			std::cerr << "Do not recognize format " << argv[i] << "." << endl;
			return 1;
		    }
		    if (second.empty()) {
			if (!first.compare("fasta") || !first.compare("phylip")) tree_file_format = "newick";
			else if (!first.compare("nexus")) tree_file_format = "nexus";
		    }
		    else if (!second.compare("nexus")) tree_file_format = "nexus";
		    else if (!second.compare("newick")) tree_file_format = "newick";
		    else {
			std::cerr << "Do not recognize format " << second << "." << endl;
			return 1;
		    }
		}
		else std::cerr << "--format require nexus or newick as additional argument" << endl;
	    }
            else if (!strcmp(argv[i],"--output")) {
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    if (!strcmp(argv[i],"nexus") || !strcmp(argv[i],"nex") || (argv[i][0] == 'x' && argv[i][1] == '\0')) print_tree = 'x';
		    else if (!strcmp(argv[i],"newick") || !strcmp(argv[i],"new") || (argv[i][0] == 'w' && argv[i][1] == '\0')) print_tree = 'w';
		    else if (!strcmp(argv[i],"no")) print_tree = 'N';
		    else {
			std::cerr << "Do not recognize format " << argv[i] << "." << endl;
			return 1;
		    }
		}
		else std::cerr << "--output require nexus(nex or x) or newick (new or w) as additional argument" << endl;
	    }
////////////////////////
            else if ((i == argc-1 && argv[i][0] != '-') || ((!strcmp(argv[i],"-f") || !strcmp(argv[i],"--file")) && i < argc-1 && argv[i+1][0] != '-' && ++i) ) {
		tree_same_as_data_when_nexus = true;
		if (data_file_name.empty() && method != 'C') data_file_name = argv[i]; // Read as data file unless datafile already given or data will be simulated
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
    if (!quiet) {
       std::cerr << "The program was called with the following command:" << endl;
       for (int i=0; i<argc; ++i) std::cerr << argv[i] << ' ';
       std::cerr << endl << endl;
    }
    #ifdef DEBUG
    cerr << "Running method: " << method << endl;
    #endif //DEBUG
    /// preparing reading from files
    if (!data_file_name.empty()) {
	#ifdef DEBUG
	cerr << data_file_name << endl;
	#endif //DEBUG
	if (!quiet) cerr << "Reading data from: " << data_file_name << endl;
        data_file.open(data_file_name.c_str(),std::ifstream::in);
        if (data_file.good())
            data_input.file_stream = &data_file;
        else {
            cerr << "Could not open file: " << data_file_name << endl;
            return 1;
        }
    }
    /*if (method == 't' && rate_mod.empty()) {
	rate_mod.push_back(1.0);
    }
    else if (method == 'A' && rate_mod.get_n_clades() < taxon_sets.size()) {
	for (unsigned int i = rate_mod.size(); i < taxon_sets.size(); ++i) {
	    rate_mod.push_back(1.0);
	}
    }*/
    /// if doing NJ do that now
    if (method == 'n') { // Neigbour Joining
	njtree tree;
	if (!quiet) cerr << "Reading distance matrix." << endl;
	tree.read_distance_matrix(*data_input.file_stream, labels);
	if (!tree.matrix_good()) {
	    cerr << "Error in distance matrix. Check distance matrix format." << endl;
	    return 1;
	}
	if (!quiet) cerr << "Read distances for " << tree.n_nodes_in_array() << " taxa." << endl;
	if (!quiet) tree.print_node_and_distance_array(cerr);
	if (!quiet) cerr << "Creating NJ tree." << endl;
	tree.build_nj_tree(quiet);
	if (print_tree == 'w') tree.print_newick(print_br_length);
	else if (print_tree == 'x') tree.print_nexus(print_br_length);
	return 0;
    }
    /// For parsimony and ML
    if (method == 'p' || method == 's' || method == 'r' || method == 'o' || method == 't' || method == 'A' || method == 'S' || method == 'C') { // ML or parsimony
	vector<character_vector> characters;
	map<char, bitset<SIZE> > alphabet;
	alphabet::set_alphabet_dna(alphabet);
	if (!alphabet_file_name.empty()) {
	    alphabet.clear();
	    if (!alphabet_file_name.compare("binary")) {
		alphabet::set_alphabet_binary(alphabet);
		if (!quiet) cerr << "Using binary alphabet." << endl;
	    }
	    else if (!alphabet_file_name.compare("dna")) {
		alphabet::set_alphabet_dna(alphabet);
		if (!quiet) cerr << "Using dna alphabet." << endl;
	    }
	    else if (!alphabet_file_name.compare("protein")) {
		alphabet::set_alphabet_amino_acid(alphabet);
		if (!quiet) cerr << "Using amino acid alphabet." << endl;
	    }
	    else {
		ifstream alphabet_file;
		alphabet_file.open(alphabet_file_name.c_str(),std::ifstream::in);
		alphabet_parser parser(alphabet_file,alphabet);
		parser.pars();
		if (!quiet) cerr << "Setting alphabet from file " << alphabet_file_name << "." << endl;
	    }
	}
	else if (!quiet) cerr << "Using default alphabet (dna)." << endl;
    	bitset<SIZE> all_possible_char;
	for (map<char,bitset<SIZE> >::const_iterator i = alphabet.begin(); i != alphabet.end(); ++i)
	    all_possible_char |= i->second;
	if (!quiet) {
	    cerr << "Alphabet include " << all_possible_char.count() << " distinct traits, given by " << alphabet.size() << " characters." << endl;
	}
	if (alphabet.empty()) {
	    cerr << "No characters in alphabet. Unable to read data." << endl;
	    return 1;
	}
	partitions regions;
	regions.add_alphabet("first",alphabet);
	regions.add_partition(0,0,"default","first");
	#ifdef DEBUG
	cerr << "Prepared partitions." << endl;
	#endif //DEBUG
	if (method != 'C') {
	    if (!data_file_format.empty()) {
		data_input.set_file_type(data_file_format.c_str());
	    }
	    else if (!data_input.set_file_type()) {
		if (!data_input.set_file_type("fasta")) cerr << "Failed to set default file type." << endl;
	    }
	    if (!quiet) cerr << "Data file format: " << data_input.get_file_type() << endl;
	    if (data_input.test_file_type("nexus")) {
		if (!data_input.move_to_next_X_block( nexus_block::DATA )) {
		    cerr << "Could not find DATA block in NEXUS file." << endl;
		    return 1;
		}
	    }
	    matrix_parser data_parser(*data_input.file_stream, characters, regions, data_input.get_file_type() );
	    data_parser.pars();
	}
	if (!quiet && !characters.empty()) {
	    cerr << "Number of taxa: " << characters.size() << "." << endl;
	    cerr << "Number of characters (first taxon): " << characters.begin()->n_char() << endl;
	}
	#ifdef DEBUG
	if (!characters.empty()) std::cerr << "Max number of characters: " << characters.begin()->max_n_char() << std::endl;
	#endif //DEBUG
        if (data_file.is_open()) data_file.close();
	if (non_as_uncertain) {
	    for (vector<character_vector>::iterator i = characters.begin(); i != characters.end(); ++i)
		i->set_empty_to(all_possible_char);
	}
	if (!tree_file_name.empty()) {
	    tree_file.open(tree_file_name.c_str(),std::ifstream::in);
	    if (tree_file.good())
		tree_input.file_stream = &tree_file;
	    else {
		cerr << "Could not open file: " << tree_file_name << endl;
		return 1;
	    }
	    if (!quiet) cerr << "Reading trees from " << tree_file_name << "." << endl;
	}
	if (tree_same_as_data_when_nexus && data_input.test_file_type("nexus") &&
	(tree_file_format.compare("nexus") || tree_file_format.empty()) &&
	tree_file_name.empty() && !data_file_name.empty()) {
    	    tree_file_name = data_file_name;
	    if (tree_file_format.empty()) tree_file_format = "nexus";
	}
    	if (!tree_file_format.empty()) {
    	    tree_input.set_file_type(tree_file_format.c_str());
	}
	else if (!tree_input.set_file_type()) {
	    if (!tree_input.set_file_type("newick")) cerr << "Failed to set default file type." << endl;
	}
	if (!quiet) cerr << "Tree file format: " << tree_input.get_file_type() << endl;
/////////////////////////////////
	if (method == 'p') {
	    if (!quiet) cerr << "Calculating parsimony scores:" << endl;
	    char nexus_command = nexus_command::NON;
    	    if (tree_input.test_file_type("nexus")) {
		if (tree_input.move_to_next_X_block( nexus_block::TREES )) {
		    nexus_command = tree_input.read_next_nexus_command();
		    if (nexus_command==nexus_command::TRANSLATE) {
			tree_input.read_translate_parameters(taxa_trans);
			nexus_command = tree_input.read_next_nexus_command();
		    }
		}
	    }
	    if (!quiet && !taxa_trans.empty()) cerr << "Read translation for " << taxa_trans.size() << " taxa." << endl;
	    unsigned int read_trees(0);
	    while (1) {
		if (tree_input.test_file_type("nexus")) {
		    if (read_trees != 0) nexus_command = tree_input.read_next_nexus_command();
		    if (!(nexus_command==nexus_command::TREE && tree_input.move_to_start_of_tree()))
			break;
		}
		tree tree;
		tree.tree_file_parser( *tree_input.file_stream, taxa_trans, false );
		if (tree.empty()) break;
		++read_trees;
		unsigned int score  = tree.fitch_parsimony( characters, print_state_on_nodelabel, calc_br_length, alphabet );
		if (print_br_length || print_state_on_nodelabel) {
		    if (print_tree == 'w') {
			cout << "[score: " << score << "] ";
			tree.print_newick(print_br_length);
		    }
		    else if (print_tree == 'x') {
			if (read_trees == 1) tree.print_nexus_tree_intro(taxa_trans);
			stringstream ss;
			ss << "tree" << read_trees;
			ss << '_' << read_trees;
			ss << " [" << score << "] ";
			string tree_name(ss.str());
			tree.print_tree_to_nexus( tree_name, print_br_length, true, taxa_trans );
		    }
		}
		else cout << score << endl;
	    }
	    if ((print_br_length || print_state_on_nodelabel) && print_tree == 'x') {
		cout << "End;" << endl;
	    }
	}
	else if (method == 's' || method == 'r') {
	    if (!quiet) cerr << "Performing stepwise addition." << endl;
	    tree tree;
	    if (method == 'r') random_shuffle(characters.begin(),characters.end());
	    tree.stepwise_addition(characters);
	    if (print_br_length) {
		if (!quiet) cerr << "Calculating branch lengths." << endl;
		tree.fitch_parsimony( characters, print_br_length );
	    }
	    if (print_tree == 'w') tree.print_newick(print_br_length);
	    else if (print_tree == 'x') tree.print_nexus(print_br_length);
	}
	else if (method == 'o' || method == 't' || method == 'A' || method == 'S' || method == 'C') { // if model based
	    if (!characters.empty() && characters.begin()->n_char() > 1) cerr << "Warning!!! Will only calculate likelihood of first character in matrix." << endl;
	    unsigned int n_states = 0;
	    //vector<double> extra_parameters;
	    if (method != 'C' && !char_freqs.empty()) tr = true;
	    if (!quiet) cerr << "Preparing model." << endl; 
	    if (method == 'o' || method == 't' || method == 'A' || method == 'S') {
		for (vector<character_vector>::iterator i=characters.begin(); i!=characters.end(); ++i)
		    if (i->highest_char_state()+1 > n_states) n_states = i->highest_char_state()+1;
	    }
	    else if (method == 'C') {
		for (map<char,bitset<SIZE> >::iterator i = alphabet.begin(); i != alphabet.end(); ++i) {
		    for (unsigned int j=n_states; j < SIZE; ++j) {
			if (i->second.test(j)) n_states = j+1;
		    }
		}
	    }
	    #ifdef DEBUG
	    cerr << "Checked number of states: " << n_states << endl;
	    #endif //DEBUG
	    sub_model model(n_states,tr);
	    #ifdef DEBUG
            cerr << "N states in model: " << model.get_n_states() << endl;
	    #endif //DEBUG
    	    if (!model_specifications.empty()) {
		// check model specifications
		for (unsigned int i=0; i < model_specifications.size(); ++i) {
		    if (!model.set_rate_spec(i, model_specifications[i])) {
			cerr << "Model miss-specification. Could not set parameter " << model_specifications[i] << " for rate " << i << "." << endl;
			cerr << "Rates need to be named in order (i.e. 1, 2, 3...)." << endl;
			return 1;
		    }
		}
	    }
    	    else {
	       	for (unsigned int i=0; model.n_posible_rates(); ++i) model.set_rate_spec(i,i);
	    }
	    #ifdef DEBUG
	    cerr << "Set model specifications." << endl;
	    #endif //DEBUG
    // Set parameter values
    	    if(!model_parameters.empty()) {
		#ifdef DEBUG
		cerr << "model parameter size: " << model_parameters.size() << ", model n rates: " << model.get_n_rates() << endl;
		#endif //DEBUG
		if (model_parameters.size() != model.get_n_rates()) {
		    cerr << "The number of parameter values (" << model_parameters.size() << ") must be as many as parameters (" << model.get_n_rates() << ")." << endl;
		    return 1;
		}
		for (unsigned int i=0; i < model_parameters.size(); ++i) {
		    #ifdef DEBUG
		    cerr << "Parameter: " << i << ", value: " << model_parameters[i] << " of " << model.get_n_rates() << " rates" << endl;
		    #endif //DEBUG
		    if (!model.set_rate(i,model_parameters[i])) {
			cerr << "All parameter values must be positive." << endl;
			return 1;
		    }
		    #ifdef DEBUG
		    cerr << "TTT Parameter: " << i << ", value: " << model_parameters[i] << " of " << model.get_n_rates() << " rates" << endl;
		    #endif //DEBUG
		}
	    }
	    else {
		for (unsigned int i=0; i < model.get_n_rates(); ++i) {
		    if (model.is_tr()) { model.set_rate(i,1.0); }
		    else model.set_rate(i,0.1);
		}
	    }
	    #ifdef DEBUG
            cerr << "Set model parameters." << endl;
            #endif //DEBUG
	    if (model.is_tr()) {
		if (char_freqs.empty()) {
		    model.set_freq_same();
		}
		else {
		    for (unsigned int i = 0; i< char_freqs.size(); ++i) {
			if(!model.set_freq(i,char_freqs[i])) {
			    cerr << "Failed to add freq " << char_freqs[i] << " for character " << i << endl;
			    if (i >= model.get_n_states()) return 1;
			    else if (i == model.get_n_states()-1) cerr << "Last char will be 1.0 minus the frequency of previous characters" << endl;
			}
		    }
		}
		#ifdef DEBUG
    		cerr << "Set state frequencies." << endl;
		#endif //DEBUG
	    }
	    // Adjust number of parameters according to rate change model and add extra parameters to vector
	    /*if (method == 't') {
		extra_parameters.push_back(cut_off);
		#ifdef DEBUG
		cerr << "N extra parameters: " << extra_parameters.size() << endl;
		#endif //DEBUG
		extra_parameters.push_back(rate_mod[0]);
		#ifdef DEBUG
		cerr << "N extra parameters: " << extra_parameters.size() << endl;
		#endif //DEBUG
	    }
	    else if (method == 'A') {
		for (unsigned int i = 0; i < taxon_sets.size(); ++i) {
		    extra_parameters.push_back(rate_mod[i]);
		    #ifdef DEBUG
		    cerr << "N extra parameters: " << extra_parameters.size() << endl;
		    #endif //DEBUG
		}
	    }*/

	    // check fixed_parameters
	    for (set<unsigned int>::iterator i=fixed_parameters.begin(); i != fixed_parameters.end(); ++i) {
		if (*i > model.get_n_rates()) {
		    cerr << "Can not fix rate parameter " << *i << ". It is out of bound (n parameters=" << model.get_n_rates() << ")." << endl;
		    return 1;
		}
	    }
	    if (method == 'S') { // if skyline fix all the cut offs
		unsigned int n = rate_mod.get_n_cut_offs ();
		set<unsigned int> temp;
		for (set<unsigned int>::iterator i=fixed_extra_parameters.begin(); i != fixed_extra_parameters.end(); ++i)  temp.insert(*i + n);
		fixed_extra_parameters = temp;
		for (; n > 0; --n) {
		    fixed_extra_parameters.insert(n-1);
		}
	    }
	    for (set<unsigned int>::iterator i=fixed_extra_parameters.begin(); i != fixed_extra_parameters.end(); ++i) {
		if (*i > rate_mod.get_n_parameters()-1) {
		    cerr << "Can not fix parameter " << *i << ". It is out of bound (n parameters=" << rate_mod.get_n_rates() << ")." << endl;
		    return 1;
		}
	    }
	    #ifdef DEBUG
            cerr << "Set fixed parameters." << endl;
            #endif //DEBUG
	    if (!quiet) {
		if (tr) cerr << "Using time reversible model" << endl;
		cerr << "Model has " << model.get_n_rates() << " rate parameters, of which " << fixed_parameters.size() << " are fixed." << endl;
		if (tr && fixed_freq) cerr << "The state frequencies are fixed." << endl;
		else if (tr) cerr << "There are " << model.get_n_states()-1 << " free state frequencies." << endl;
	    }
///////////////
	    char nexus_command = nexus_command::NON;
    	    if (tree_input.test_file_type("nexus")) {
		if (tree_input.move_to_next_X_block( nexus_block::TREES )) {
		    nexus_command = tree_input.read_next_nexus_command();
		    if (nexus_command==nexus_command::TRANSLATE) {
			tree_input.read_translate_parameters(taxa_trans);
			nexus_command = tree_input.read_next_nexus_command();
		    }
		}
	    }
	    if (!quiet && !taxa_trans.empty()) cerr << "Read translation for " << taxa_trans.size() << " taxa." << endl;
	    unsigned int read_trees(0);
	    if (!quiet && method != 'C') cerr << "Estimating likelihood." << endl;
	    #ifdef DEBUG
	    cerr << "Starting to read trees." << endl;
	    #endif //DEBUG
	    while (1) {
		stringstream model_out;
		if (tree_input.test_file_type("nexus")) {
		    if (read_trees != 0) nexus_command = tree_input.read_next_nexus_command();
		    if (!(nexus_command==nexus_command::TREE && tree_input.move_to_start_of_tree()))
			break;
		}
		#ifdef DEBUG
		cerr << "Atempt to read tree." << endl;
		#endif //DEBUG
		simpleML tree;
		tree.tree_file_parser( *tree_input.file_stream, taxa_trans, false );
		if (tree.empty()) break;
		++read_trees;
		if (!quiet) cerr << "Read tree " << read_trees << endl;
		if (!taxon_sets.empty()) tree.add_taxon_sets(taxon_sets);
		if (tree.shortest_branch() < 0.0) cerr << "Warning!!! Tree " << read_trees << " contains negative branches." << endl;
		tree.init(n_states); // need to fix this
		for (vector<character_vector>::iterator i=characters.begin(); i!=characters.end(); ++i)
		    tree.set_char(i->get_taxon(),i->get_character(0));
////////////////////////////////////////////
		if (print_tree == 'x') model_out << '[' << endl;
		if (method == 'C') {
		    #ifdef DEBUG
		    cerr << "simulating data" << endl;
		    for (int i=0; i < char_freqs.size(); ++i) cerr << char_freqs[i] << endl;
		    #endif //DEBUG
		    characters = tree.simulate_chars(char_freqs, n_char, model);

		}
		else if (optimize_param) {
		    #ifdef NLOPT
		    if (!quiet) cerr << "Optimizing parameters" << endl;
		    vector<double> lower_bounds;
		    vector<double> upper_bounds;
		    vector<double> variable_values; // starting values for variables to optimize
		    unsigned int start_freq(0);
		    unsigned int start_extras(0);
		    unsigned int n_parameters(0);
		    // add substitution rates to variables to optimize
		    for (unsigned int i=0; i < model.get_n_rates(); ++i) {
			if (fixed_parameters.find(i) == fixed_parameters.end()) {
			    variable_values.push_back(model.get_rate(i));
			    lower_bounds.push_back(0.0);
			    upper_bounds.push_back(DBL_MAX);
			    ++n_parameters;
			}
		    }
		    // add frequencies to model parameters to optimize
		    if (tr && !fixed_freq) {
			start_freq = n_parameters;
			for (unsigned int i=0; i < model.get_n_states()-1; ++i) {
    			    lower_bounds.push_back(0.0);
			    upper_bounds.push_back(1.0);
			    variable_values.push_back(model.get_freq(i));
			    ++n_parameters;
	    		}
		    }
		    // add time and rate modifiers to parameters to optimize
		    if (rate_mod.get_n_rates()) {
			start_extras = n_parameters;
			for (unsigned int i=0; i < rate_mod.get_n_parameters(); ++i) {
			    if (fixed_extra_parameters.find(i) == fixed_extra_parameters.end()) {
				lower_bounds.push_back(0.0);
				if ((method == 't' || method == 'S') && i < rate_mod.get_n_cut_offs()) {
				    upper_bounds.push_back(tree.longest_to_tip());
				    variable_values.push_back(rate_mod.get_cut_off(0));
				    ++n_parameters;
				}
				else {
				    upper_bounds.push_back(DBL_MAX);
				    if (method == 't' || method == 'S') variable_values.push_back(rate_mod.get_rate_in_time(i-rate_mod.get_n_cut_offs()));
				    else if (method == 'A') variable_values.push_back(rate_mod.get_rate_clade(i));
				    ++n_parameters;
				}
			    }
	    		}
		    }
		    if (n_parameters < 1) {
			cerr << "No parameters to optimize. Try running with --no_optim." << endl;
			return 1;
		    }
		    //// Prep model
		    nlopt::opt maximize(nlopt::LN_NELDERMEAD, n_parameters);
		    tree_modelspec_struct data = {&tree, &model, &fixed_parameters, fixed_freq, start_freq, &rate_mod, start_extras, &fixed_extra_parameters};
		    maximize.set_lower_bounds(lower_bounds);
		    maximize.set_upper_bounds(upper_bounds);
		    ///////////////////
		    double LogL;
		    cerr << method << endl;
		    if (method == 'o') {
			if (variable_values.size() != maximize.get_dimension()) {
			    cerr << "Mismatch in number of starting values ("<< variable_values.size() << ") and number of parameters to optimize (" << maximize.get_dimension() << "). Quiting." << endl;
			    return 1;
			}
			maximize.set_max_objective(opt_function,&data);
			model_out << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			nlopt::result result = maximize.optimize(variable_values,LogL);
			change_non_fixed(variable_values,&data);
			if (result < 0) {
			    cerr << "Failure when optimizing!!!" << endl;
			    return 1;
			}
		    }
		    else if (method == 't' || method == 'S') {
			if (!quiet) cerr << "Estimating rates through time" << endl;
			maximize.set_max_objective(opt_rate_in_time,&data);
			if (variable_values.size() != maximize.get_dimension()) {
			    cerr << "Mismatch in number of starting values ("<< variable_values.size() << ") and number of parameters to optimize (" << maximize.get_dimension() << "). Quitting." << endl;
			    return 1;
			}
			#ifdef DEBUG
			cerr << "Fixed extra parameters" << endl;
			for (set<unsigned int>::const_iterator i = fixed_extra_parameters.begin(); i != fixed_extra_parameters.end(); ++i) cerr << *i << ' ';
			cerr << endl << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			#endif //DEBUG
			model_out << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			nlopt::result result = maximize.optimize(variable_values,LogL);
			#ifdef DEBUG
			cerr << "Optimized!!!" << endl;
			#endif //DEBUG
			if (result < 0) {
			    cerr << "Failure when optimizing!!!" << endl;
			    return 1;
			}
			change_non_fixed(variable_values,&data);
			for (unsigned int i=0; i < rate_mod.get_n_cut_offs(); ++i) {
			    model_out << "Branch length multiplied by: " << rate_mod.get_rate_in_time(i) << " until " <<  data.rate_mod->get_cut_off(i) << endl;
			}
			/*model_out << "Time from root to rate change: " << data.rate_mod->get_cut_off(0) << endl;
			model_out << "Rate multiplier: " << data.rate_mod->get_rate_in_time(0) << endl;*/
		    }
		    else if (method == 'A') {
			if (variable_values.size() != maximize.get_dimension()) {
			    cerr << "Mismatch in number of starting values ("<< variable_values.size() << ") and number of parameters to optimize (" << maximize.get_dimension() << "). Quitting." << endl;
			    return 1;
			}
			#ifdef DEBUG
			cerr << "Variable values: ";
			for (vector<double>::const_iterator i=variable_values.begin(); i != variable_values.end(); ++i) cerr << *i << ' ';
			cerr << endl;
			#endif //DEBUG
			maximize.set_max_objective(opt_rate_for_clades_function,&data);
			model_out << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			nlopt::result result = maximize.optimize(variable_values,LogL);
			if (result < 0) {
			    cerr << "Failure when optimizing!!!" << endl;
			    return 1;
			}
			#ifdef DEBUG
			cerr << "Variable values 2: ";
			for (vector<double>::const_iterator i=variable_values.begin(); i != variable_values.end(); ++i) cerr << *i << ' ';
			cerr << endl;
			#endif //DEBUG
			change_non_fixed(variable_values,&data);
			/// need to fix if any rates are fixed
			for (unsigned int i=0; i < taxon_sets.size(); ++i) {
			    //model_out << "Rate mod for \"" << taxon_sets[i] << "\": " << variable_values[data.extra_start+i] << endl;
			    model_out << "Rate mod for (" << i << ") \"" << taxon_sets[i] << "\": " << rate_mod.get_rate_clade(i) << endl;
			}
		    }
		    else {
			cerr << "Unrecognized method. Nothing to do." << endl;
			return 1;
		    }
		    #if DEBUG
		    model_out << "Parameter values:" << endl;
		    for (unsigned int i= 0; i < variable_values.size(); ++i) {
			model_out << i << ": " << variable_values[i] << ' ';
		    }
		    model_out << endl;
		    #endif /*DEBUG*/
		    if (tr) {
			model_out << "Rates:" << endl;
			for (unsigned int i=0; i < n_states; ++i) {
			    bitset<SIZE> I;
			    I.set(i);
			    for (unsigned int j=i+1; j < n_states; ++j) {
				bitset<SIZE> J;
				J.set(j);
				if (i != 0 || j != i+1) model_out << "; ";
				model_out << alphabet::translate_bitset(I, alphabet) << "<->" << alphabet::translate_bitset(J, alphabet) << ": ";
				model_out << model.get_rate(i,j);
			    }
			}
			model_out << endl;
			model_out << "Frequencies:" << endl;
			for (unsigned int i=0; i < n_states; ++i) {
			    if (i != 0) model_out << "; ";
			    bitset<SIZE> temp;
			    temp.set(i);
			    model_out << alphabet::translate_bitset(temp, alphabet) << ": " << model.get_freq(i);

			}
			model_out << endl;
		    }
		    model_out << "Q matrix:" << endl;
		    model.print_Q_matrix( model_out );
		    model_out << endl;
		    model_out << "Log likelihood: " << LogL << endl;
/////////////////////
		    #else
		    cerr << "Compiled without NLOPT can not do optimization" << endl;
		    #endif //NLOPT
		}
		/////////////////////////////
		else {
		    if (method == 'o') {
			model_out << "Log likelihood= " << tree.calculate_log_likelihood(model) << endl;
		    }
		    else if (method == 't') {
			model_out << "Log likelihood= " << tree.calculate_likelihood_rate_change_in_time(rate_mod,model) << endl;
			model_out << "Time from root: " << rate_mod.get_cut_off(0) << " Rate modifier: " << rate_mod.get_rate_in_time(0) << endl;
		    }
		    else if (method == 'A') {
			model_out << "Log likelihood= " << tree.calculate_likelihood_rate_change_at_nodes ( rate_mod, model ) << endl;
			for (unsigned int i=0; i < taxon_sets.size(); ++i) {
			    model_out << "Rate mod for \"" << taxon_sets[i] << "\": " <<  rate_mod.get_rate_clade(i) << endl;
			}
		    }
		    else {
			cerr << "Unrecognized method. Nothing to do." << endl;
			return 1;
		    }
		    if (tr) {
			model_out << "Rates:" << endl;
			for (unsigned int i=0; i < n_states; ++i) {
			    bitset<SIZE> I;
			    I.set(i);
			    for (unsigned int j=i+1; j < n_states; ++j) {
				bitset<SIZE> J;
				J.set(j);
				if (i != 0 || j != i+1) model_out << "; ";
				model_out << alphabet::translate_bitset(I, alphabet) << "<->" << alphabet::translate_bitset(J, alphabet) << ": ";
				model_out << model.get_rate(i,j);
			    }
			}
			model_out << endl;
			model_out << "Frequencies:" << endl;
			for (unsigned int i=0; i < n_states; ++i) {
			    if (i != 0) model_out << "; ";
			    bitset<SIZE> temp;
			    temp.set(i);
			    model_out << alphabet::translate_bitset(temp, alphabet) << ": " << model.get_freq(i);

			}
			model_out << endl;
		    }
		    else {
			model_out << "Parameters:" << endl;
			model.print_Q_matrix( model_out );
		    }
		    model_out << endl;
		}
		if (method == 'C') {
		    for (vector<character_vector>::iterator i = characters.begin(); i != characters.end(); ++i) {
			cout << '>' << i->get_taxon() << endl;
			for (unsigned int j =0; j < i->n_char(); ++j) {
			    cout << alphabet::translate_bitset(i->get_character(j),alphabet);
			}
			cout << endl;
		    }
		}
		else if (print_tree == 'w') {
		    // cout << "Tree:" << endl;
		    if (print_state_on_nodelabel) tree.draw_normalized_likelihood_on_nodes();
		    tree.print_newick();
		}
		else if (print_tree == 'x') {
		    model_out << ']' << std::endl;
		    if (read_trees == 1) tree.print_nexus_tree_intro(taxa_trans);
		    if (print_state_on_nodelabel) tree.draw_normalized_likelihood_on_nodes();
		    stringstream ss;
		    ss << "tree" << read_trees;
		    ss << '_' << read_trees;
		    string tree_name(ss.str());
		    tree.print_tree_to_nexus( tree_name, print_br_length, true, taxa_trans );
		    //tree.print_tree_to_nexus();
		}
		else {}
		cout << model_out.str();
/////////////////////////////
	    }
	    if (print_tree == 'x') {
		cout << "End;" << endl;
	    }
	}
    }
}

void help () {
    cout << "Treeator " << VERSION <<  " is a command line program to construct trees. The program take" << endl;
    cout << "either a left triangular similarity matrix (neighbour joining) or a data matrix" << endl;
    cout << "of fasta, nexus, or relaxed phylip format (not interleaved; parsimony/maximum" << endl;
    cout << "likelihood) as input through standard in/ last argument/ as given below. For" << endl;
    cout << "character data an alphabet is also needed, default is binary (0 1 -)." << endl;
    cout << "(c) Martin Ryberg " << YEAR << "." << endl << endl;
    cout << "Usage:" << endl << "treeator [arguments] data_file.txt" << endl << endl;
    cout << "Arguments:" << endl;
    cout << "--alphabet_file / -a [file/type] give file with character alphabet, or dna," << endl;
    cout << "                                 protein, or binary for dna, amino acid," << endl;
    cout << "                                 respectively binary (0 1) alphabets (default:" << endl;
    cout << "                                 dna)." << endl;
    cout << "--clade_rates / -V                give clades by giving a sets of taxa that have" << endl;
    cout << "                                 the base of the clade as most recent common" << endl;
    cout << "                                 ancestor. Different sets should be separated by" << endl;
    cout << "                                 semicolon (;), taxa in a set should be separated" << endl;
    cout << "                                 by comma. The rate of the set should be given" << endl;
    cout << "                                 before the set and separated from it by a colon" << endl;
    cout << "                                 (:). If a set is not given a rate it will get" << endl;
    cout << "                                 the same rate as the previous set, and if" << endl;
    cout << "                                 optimizing the rate parameters they will be" << endl;
    cout << "                                 optimized as one. E.g. -A \"0.5:Taxon1,Taxon2;" << endl;
    cout << "                                 Taxon3,Taxon4\"" << endl;
    cout << "--data_file / -d [file]          give the data file."<< endl;
    cout << "--do_not_calc_br_length / -c     do not save branch length when doing parsimony" << endl;
    cout << "                                 reconstruction. This means that the branch" << endl;
    cout << "                                 lengths of the input tree are kept. If you do" << endl;
    cout << "                                 not want branch length in the output tree use" << endl;
    cout << "                                 --no_branch_length/-0." << endl;
    cout << "--fixed / -e [number/s]          give parameter to fix. First parameter is" << endl;
    cout << "                                 indexed 0. Several parameters can be given in a" << endl;
    cout << "                                 comma separated string, e.g. -e 0,2,3." << endl;
    cout << "--fix_freq                       set state frequences to be equal" << endl;
    cout << "--fixed_extras                   give if rate change parameters should be fixed" << endl;
    cout << "                                 (under -V, -S, and -U). Under -U the cut off in" << endl;
    cout << "                                 time for rate change will be parameter 0 and the" << endl;
    cout << "                                 rate multiplier parameter 1. e.g. --fixed_extras" << endl;
    cout << "                                 0" << endl;
    cout << "--file / -f [file]               give data file name, or if data file name" << endl;
    cout << "                                 already given, then tree file name. If nexus" << endl;
    cout << "                                 format and no tree file name is given, tree is" << endl;
    cout << "                                 assumed to be in same file as data." << endl;
    cout << "--format [format]                give the format of the input files. For" << endl;
    cout << "                                 character file fasta, phylip and nexus are the" << endl;
    cout << "                                 options. For the tree file the options are" << endl;
    cout << "                                 newick and nexus. Give the character file" << endl;
    cout << "                                 format first, and the tree file format after a" << endl;
    cout << "                                 comma, e.g. --format phylip,newick. Fasta is" << endl;
    cout << "                                 the default character file format, and newick" << endl;
    cout << "                                 is the default tree file format (unless" << endl;
    cout << "                                 character file is set to nexus, then nexus is" << endl;
    cout << "                                 also default for tree file, e.g. --format" << endl;
    cout << "                                 nexus)." << endl;
    cout << "--frequencies / -F               give state frequencies as comma separated" << endl;
    cout << "                                 string, e.g. -F 0.25,0.25,0.25,0.25" << endl;
    cout << "--get_state_at_nodes             will give the states at internal nodes as" << endl;
    cout << "                                 comments (readable in FigTree)." << endl;
    cout << "--help / -h                      print this help." << endl;
    cout << "--likelihood / -l                calculate likelihood for data given tree." << endl;
    cout << "--model / -m [number/s]          give the model by numbering the rate parameters" << endl;
    cout << "                                 for different transition, e.g. -m 0,1,0,2,1,2." << endl;
    cout << "                                 The order is by row, i.e. from parameter 0 to" << endl;
    cout << "                                 parameter 1 first then, 0 to 2, and so on to" << endl;
    cout << "                                 all other parameters, then 1 to 0, and so" << endl;
    cout << "                                 forth." << endl;
    cout << "--neighbour_joining / -n         compute neighbour joining tree for given data." << endl;
    cout << "                                 The data should be  a left triangular" << endl;
    cout << "                                 similarity matrix." << endl;
    cout << "--no_branch_length / -0          do not print branch lengths and do not" << endl;
    cout << "                                 calculate branch lengths for parsimony trees." << endl;
    cout << "--no_label / -L                  will tell treeator that there are no taxon" << endl;
    cout << "                                 labels in the similarity matrix." << endl;
    #ifdef NLOPT
    cout << "--no_optim / -N                  calculate likelihood for given parameters. No" << endl;
    cout << "                                 optimization." << endl;
    #endif //NLOPT
    cout << "--output [newick/nexus]          give tree format for output, nexus (nex or x" << endl;
    cout << "                                 for short) or newick (new or w for short), e.g" << endl;
    cout << "                                 --output x. (default w)." << endl; 
    cout << "--parameters / -P [values]       give corresponding parameter values for" << endl;
    cout << "                                 parameters.";
    #ifdef NLOPT
    cout << " If optimizing these will be" << endl;
    cout << "                                 starting values, e.g.";
    #else
    cout << " E.g.";
    #endif //NLOPT
    cout << " -P 0.1,0.01,0.05.";
    cout << endl;
    cout << "--parsimony / -p                 calculate parsimony score for given tree and" << endl;
    cout << "                                 data." << endl;
    cout << "--random / -r                    do stepwise addition in random order." << endl;
    /*cout << "--rate_mod / -R [value]          give modifier for rate compared to rate at root" << endl;
    cout << "                                 e.g. -r 0.5. Default: 1.0." << endl;*/
    cout << "--rate_change_at_time / -U       set a point in time, from the root, during" << endl;
    cout << "                                 which the rate will be multiplied by a number." << endl;
    cout << "                                 The rate multiplier should be given first" << endl;
    cout << "                                 folowed by a colon (:) and then the time, e.g." << endl;
    cout << "                                 -U 0.5:50" << endl;
    cout << "--simulate                       simulate data on given tree. Number of" << endl;
    cout << "                                 characters as possible extra argument (default" << endl;
    cout << "                                 1). Use model given by -m and parameters givens" << endl;
    cout << "                                 by -P and probability of each character at roots" << endl;
    cout << "                                 by -F, e.g. --simulate 10." << endl;
    /*cout << "--time / -T [value]              give branch length distance from root where" << endl;
    cout << "                                 change in rate occur, e.g. -T 10. Default: 0." << endl;*/
    cout << "--time_reversible / --tr         set the model to be time reversible. If no" << endl; 
    cout << "                                 frequencies (-F) are given. State frequencies" << endl;
    cout << "                                 are assumed to be equal" << endl; 
    cout << "--tree_file / -t [file]          give tree file name." << endl;
    cout << "--skyline_rates / -S             add rate multipliers up until time points. Give" << endl;
    cout << "                                 a vector of rates/rate staring values separateds" << endl;
    cout << "                                 from the time breaking points by a colon (:)." << endl;
    cout << "                                 each rate that should be estimated separatly" << endl;
    cout << "                                 should be given separatly. Time points are" << endl;
    cout << "                                 fixed. E.g. 0.5:50,110,2.0:75. Each rate will" << endl;
    cout << "                                 be from previus time point up until the next," << endl;
    cout << "                                 times are given from the root." << endl;
    cout << "--step_wise / -s                 do parsimony stepwise addition." << endl;
    cout << "--verbose / -v                   get additional output." << endl;
    cout << endl;
}

double opt_function(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    if (!change_non_fixed(x,static_cast<tree_modelspec_struct*>(data))) return 0.0;
    double LogLH = static_cast<tree_modelspec_struct*>(data)->tree->calculate_log_likelihood(*(static_cast<tree_modelspec_struct*>(data)->model));
    #ifdef DEBUG
        std::cerr << "Likelihood: " << LogLH << endl;
    #endif //DEBUG
    return LogLH;
}

double opt_rate_in_time(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    if (!change_non_fixed(x,static_cast<tree_modelspec_struct*>(data))) return 0.0;
    if (static_cast<tree_modelspec_struct*>(data)->rate_mod->get_n_time_rates() < 1) return 0.0;
    double LogLH =  static_cast<tree_modelspec_struct*>(data)->tree->calculate_likelihood_rate_change_in_time(*(static_cast<tree_modelspec_struct*>(data)->rate_mod),*(static_cast<tree_modelspec_struct*>(data)->model));
    return LogLH;
}

double opt_rate_for_clades_function(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    if (!change_non_fixed(x,static_cast<tree_modelspec_struct*>(data))) return 0.0;
    #ifdef DEBUG
    cerr << "Rates:";
    for (unsigned int i=0; i < static_cast<tree_modelspec_struct*>(data)->rate_mod->get_n_clade_rates(); ++i) { cerr << " " << static_cast<tree_modelspec_struct*>(data)->rate_mod->get_rate_clade(i); }
    cerr << endl;
    #endif //DEBUG
    // create a new rate_mod vector based on specs
    double LogLH = static_cast<tree_modelspec_struct*>(data)->tree->calculate_likelihood_rate_change_at_nodes(*(static_cast<tree_modelspec_struct*>(data)->rate_mod),*(static_cast<tree_modelspec_struct*>(data)->model));
    return LogLH;

}


bool change_non_fixed(const std::vector<double> &x, tree_modelspec_struct* data) {
    unsigned int j(0);
    // Set substitution model parameters
    for (unsigned int i=0; i < data->model->get_n_rates(); ++i) {
	if (data->fixed->find(i) == data->fixed->end()) {
	    if (!data->model->set_rate(i,x[j])) return false;
	    #ifdef DEBUG
	    cerr << "Set rate " << i << " to " << x[j] << endl;
	    #endif //DEBUG
	    ++j;
	}
    }
    // if time reversable model, and state frequences not fixed set state frequencies
    if (data->model->is_tr() && !data->fixed_freq) {
	for (unsigned int i=0; i < data->model->get_n_states()-1; ++i) {
	    if (!data->model->set_freq(i,x[data->freq_start+i])) return false;
	}
    }
    // set rate model parameters
    j = 0;
    for (unsigned int i = 0; i < data->rate_mod->get_n_parameters(); ++i) {
	if (data->fixed_extra->find(i) == data->fixed_extra->end()) {
	    #ifdef DEBUG
	    cerr << "Set rate model parameter " << i << " to " << x[data->extra_start+j] << endl;
	    #endif //DEBUG
	    data->rate_mod->set_parameter(i,x[data->extra_start+j]);
	    ++j;
	}
    }
    return true;
}

