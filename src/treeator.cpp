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

using namespace std;

void help ();

double opt_function(const std::vector<double> &x, std::vector<double> &grad, void* data);
double opt_rate_in_time(const std::vector<double> &x, std::vector<double> &grad, void* data);
double opt_rate_for_clades_function(const std::vector<double> &x, std::vector<double> &grad, void* data);
void change_non_fixed(const std::vector<double> &x, std::vector<double>* values, const std::set<unsigned int>* fixed);

struct tree_modelspec_struct {
    simpleML* tree;
    const vector<unsigned int>* specifications; // which parameter represents which rate i.e. 0,0,1,1,2,2 represents same rates going from 0->1 as 0->2, and 1->0 as 1->2, and so on.
    vector<double>* values;
    const set<unsigned int>* fixed;
    bool tr;
};

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
    double cut_off = 0.0;
    vector <double> rate_mod;
    set<unsigned int> fixed_parameters;
    vector<double> model_parameters;
    vector<double>char_freqs;
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
		    while (!numbers.empty()) {
                        char_freqs.push_back(atof(numbers.back().c_str()));
			numbers.pop_back();
		    }
		}
                else {
                    cerr << "-F/--frequencies require a comma separated real number string as next argument.";
                    return 1;
                }
		tr = true;
		if (method != 'A' && method != 't' && method != 'C') method = 'o';
	    }
            else if (!strcmp(argv[i],"-N") || !strcmp(argv[i],"--no_optim")) optimize_param = false;
            else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--likelihood")) {
		    if (method != 't' && method != 'A') method = 'o';
            }
            else if (!strcmp(argv[i],"-R") || !strcmp(argv[i],"--rate_mod")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    vector<string> temp;
		    argv_parser::pars_sub_args(argv[++i], ';', temp);
		    for (vector<string>::iterator t=temp.begin(); t != temp.end(); ++t) {
			rate_mod.push_back(atof(t->c_str()));
		    }
		}
                //    rate_mod = atof(argv[++i]);
                else {
                    std::cerr << "-R/--rate_mod require a number as next argument." << endl;
                    return 1;
                }
		if (method != 'A') method='t';
            }
            else if (!strcmp(argv[i],"-T") || !strcmp(argv[i],"--time")) {
                if ( i < argc-1 && argv[i+1][0] != '-' )
                    cut_off = atof(argv[++i]);
                else {
                    std::cerr << "-T/--time require a number as next argument." << endl;
                    return 1;
                }
		method='t';
            }
	    else if (!strcmp(argv[i],"-A") || !strcmp(argv[i],"--taxon_sets")) {
		if ( i < argc-1 && argv[i+1][0] != '-' )
		    argv_parser::pars_sub_args(argv[++i],';',taxon_sets);
		else {
		    cerr << "-A / --taxon_sets require at least one set of taxa given as a comma separated string as next argument." << endl;
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
    if (method == 't' && rate_mod.empty()) {
	rate_mod.push_back(1.0);
    }
    else if (method == 'A' && rate_mod.size() < taxon_sets.size()) {
	for (int i = rate_mod.size(); i < taxon_sets.size(); ++i) {
	    rate_mod.push_back(1.0);
	}
    }
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
	//tree.print_node_and_distance_array();
	if (!quiet) cerr << "Creating NJ tree." << endl;
	tree.build_nj_tree(quiet);
	if (print_tree == 'w') tree.print_newick(print_br_length);
	else if (print_tree == 'x') tree.print_nexus(print_br_length);
	return 0;
    }
    /// For parsimony and ML
    if (method == 'p' || method == 's' || method == 'r' || method == 'o' || method == 't' || method == 'A' || method == 'C') { // ML or parsimony
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
	else if (method == 'o' || method == 't' || method == 'A' || method == 'C') {
	    if (!characters.empty() && characters.begin()->n_char() > 1) cerr << "Warning!!! Will only calculate likelihood of first character in matrix." << endl;
	    unsigned int n_states = 0;
	    unsigned int n_parameters = 0;
	    if (!quiet) cerr << "Preparing model." << endl; 
	    if (method == 'o' || method == 't' || method == 'A') {
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
    	    if (!model_specifications.empty()) {
		// check model specifications
		vector<unsigned int> temp = model_specifications;
		sort(temp.begin(),temp.end());
		for (vector<unsigned int>::const_iterator i=temp.begin(); i!=temp.end(); ++i) {
		    if (i == temp.begin()) {
			if (*i != 0) {
			    cerr << "Lowest parameter index in model (given by -m / --model) should be 0 not " << *i << "." << endl;
			    return 1;
			}
		    }
		    else {
			if (*i - *(i-1) > 1) { 
			    cerr << "The parameter indexes should be a continuous series of integers. Missing index " << *(i-1)+1 << "." << endl;
			    return 1;
			}
		    }
		}
		set<unsigned int> parameters;
		for (vector<unsigned int>::const_iterator i=model_specifications.begin(); i!=model_specifications.end(); ++i) {
		    parameters.insert(*i);
		}
		n_parameters=parameters.size();
	    }
    	    else { 
		if (tr) {
		    for (unsigned int i=0; i < (((n_states*n_states)-n_states)/2)+n_states-1; ++i) model_specifications.push_back(i);
		}
		else {
		    for (unsigned int i=0; i < (n_states*n_states)-n_states; ++i) model_specifications.push_back(i);
		}
		n_parameters=model_specifications.size();
	    }
    // Set parameter values
    	    if(!model_parameters.empty()) {
		if (model_parameters.size() != n_parameters) {
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
		if (tr) {
		    for (unsigned int i=0; i < n_states-1; ++i) model_parameters.push_back(1.0/n_states);
		}
		else {
		    for (unsigned int i=0; i < n_parameters; ++i) model_parameters.push_back(0.1);
		}
	    }
	    // Adjust number of parameters according to rate change model and add extra parameters to vector
	    #ifdef DEBUG
	    cerr << "N parameters: " << n_parameters << ", method=" << method << endl;
	    #endif //DEBUG
	    if (method == 't') {
		n_parameters+=2;
		model_parameters.push_back(cut_off);
		model_parameters.push_back(rate_mod[0]);
	    }
	    else if (method == 'A') {
		n_parameters += taxon_sets.size();
		for (int i = 0; i < taxon_sets.size(); ++i) {
		    model_parameters.push_back(rate_mod[i]);
		}
	    }
	    // check fixed_parameters
	    #ifdef DEBUG
	    cerr << "N parameters: " << n_parameters << ", method=" << method << endl;
	    #endif //DEBUG
	    for (set<unsigned int>::iterator i=fixed_parameters.begin(); i != fixed_parameters.end(); ++i) {
		if (*i > n_parameters-1) {
		    cerr << "Can not fix parameter " << *i << ". It is out of bound (n parameters=" << n_parameters << ")." << endl;
		    return 1;
		}
	    }
	    if (!quiet) 
		cerr << "Model has " << n_parameters << " parameters, of which " << fixed_parameters.size() << " are fixed." << endl;
	    n_parameters-=fixed_parameters.size();
	    if (n_parameters==0)
		optimize_param=false;


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
		    tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
		    characters = tree.simulate_chars(char_freqs, n_char);

		}
		else if (optimize_param) {
		    #ifdef NLOPT
		    nlopt::opt maximize(nlopt::LN_NELDERMEAD, n_parameters);
		    tree_modelspec_struct data = {&tree, &model_specifications, &model_parameters, &fixed_parameters,tr};
		    vector<double> variable_values;
		    for (unsigned int i=0; i<model_parameters.size(); ++i) {
			if (fixed_parameters.find(i)==fixed_parameters.end()) variable_values.push_back(model_parameters[i]);
		    }
		    vector<double> lower_bounds;
		    for (unsigned int i=0; i < n_parameters; ++i) lower_bounds.push_back(0.0);
		    maximize.set_lower_bounds(lower_bounds);
		    double LogL;
		    //if (print_tree == 'x') model_out << '[';
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
		    else if (method == 'A') {
			maximize.set_max_objective(opt_rate_for_clades_function,&data);
			model_out << "Number of parameters to optimize: " << maximize.get_dimension() << endl;
			nlopt::result result = maximize.optimize(variable_values,LogL);
			if (result < 0) {
			    cerr << "Failure when optimizing!!!" << endl;
			    return 1;
			}
			for (int i=0; i < taxon_sets.size(); ++i) {
			    model_out << "Rate mod for \"" << taxon_sets[i] << "\": " << model_parameters[model_parameters.size()-taxon_sets.size()+i] << endl;
			}
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
		    #else
		    cerr << "Compiled without NLOPT can not do optimization" << endl;
		    #endif //NLOPT
		}
		/////////////////////////////
		else {
		    if (method == 'o') {
			tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
			model_out << "Log likelihood= " << tree.calculate_log_likelihood() << endl;
		    }
		    else if (method == 't') {
			tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
			model_out << "Log likelihood= " << tree.calculate_likelihood_rate_change_in_time(cut_off,rate_mod[0]) << endl;
			model_out << "Time from root: " << cut_off << " Rate modifier: " << rate_mod[0] << endl;
		    }
		    else if (method == 'A') {
			tree.set_Q_matrix(&model_specifications[0],&model_parameters[0]);
			model_out << "Log likelihood= " << tree.calculate_likelihood_rate_change_at_nodes ( rate_mod.data() ) << endl;
			for (int i=0; i < taxon_sets.size(); ++i) {
			    model_out << "Rate mod for \"" << taxon_sets[i] << "\": " <<  rate_mod[i] << endl;
			}
		    }
		    else {
			cerr << "Unrecognized method. Nothing to do." << endl;
			return 1;
		    }
		    tree.print_Q_matrix( model_out );
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
    cout << "--data_file / -d [file]          give the data file."<< endl;
    cout << "--do_not_calc_br_length / -c     do not save branch length when doing parsimony" << endl;
    cout << "                                 reconstruction. This means that the branch" << endl;
    cout << "                                 lengths of the input tree are kept. If you do" << endl;
    cout << "                                 not want branch length in the output tree use" << endl;
    cout << "                                 --no_branch_length/-0." << endl;
    cout << "--fixed / -e [number/s]          give parameter to fix. First parameter is" << endl;
    cout << "                                 indexed 0. Several parameters can be given in a" << endl;
    cout << "                                 comma separated string, e.g. -e 0,2,3." << endl;
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
    cout << "--get_state_at_nodes             will give the states at internal nodes as" << endl;
    cout << "                                 comments (readable in FigTree)." << endl;
    cout << "--help / -h                      print this help." << endl;
    cout << "--likelihood / -l                calculate likelihood for data given tree." << endl;
    /*cout << "--likelihood / -l [const/time]   calculate likelihood for data given tree." << endl;
    cout << "                                 Either with constant rate through time (const)" << endl;
    cout << "                                 or with rate changing (multiplied by a" << endl;
    cout << "                                 variable) at a certain time point (time)." << endl;*/
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
    /*cout << "--rate_mod / -R [value]        give modifier for rate compared to rate at root" << endl;
    cout << "                                 e.g. -r 0.5. Default: 1.0." << endl;*/
    cout << "--random / -r                    do stepwise addition in random order." << endl;
    cout << "--simulate                       simulate data on given tree. Number of" << endl;
    cout << "                                 characters as possible extra argument (default" << endl;
    cout << "                                 1). Use model given by -m and parameters givens" << endl;
    cout << "                                 by -P and probability of each character at roots" << endl;
    cout << "                                 by -F, e.g. --simulate 10." << endl;
    /*cout << "--time / -T [value]            give branch length distance from root where" << endl;
    cout << "                                 change in rate occur, e.g. -T 10. Default: 0." << endl;*/
    cout << "--tree_file / -t [file]          give tree file name." << endl;
    /*cout << "--taxon_sets / -A                give sets of taxa. Different sets should be" << endl;
    cout << "                                 separated by semicolon (;), taxa in set should" << endl;
    cout << "                                 be separated by comma, e.g. -A \"Taxon1,Taxon2;" << endl;
    cout << "                                 Taxon3,Taxon4\"." << endl;*/
    cout << "--step_wise / -s                 do parsimony stepwise addition." << endl;
    cout << "--verbose / -v                   get additional output." << endl;
    cout << endl;
}

double opt_function(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    change_non_fixed(x,static_cast<tree_modelspec_struct*>(data)->values,static_cast<tree_modelspec_struct*>(data)->fixed);
    if (static_cast<tree_modelspec_struct*>(data)->tr) static_cast<tree_modelspec_struct*>(data)->tree->set_Q_matrix(&(static_cast<tree_modelspec_struct*>(data)->specifications->front()),&(static_cast<tree_modelspec_struct*>(data)->values->at(0)),&(static_cast<tree_modelspec_struct*>(data)->values->at(static_cast<tree_modelspec_struct*>(data)->tree->get_n_states()-1))); 
    else static_cast<tree_modelspec_struct*>(data)->tree->set_Q_matrix(&(static_cast<tree_modelspec_struct*>(data)->specifications->front()),&(static_cast<tree_modelspec_struct*>(data)->values->at(0)));
    double LogLH = static_cast<tree_modelspec_struct*>(data)->tree->calculate_log_likelihood();
    //double LogLH = static_cast<tree_modelspec_struct*>(data)->tree->calculate_log_likelihood(&(static_cast<tree_modelspec_struct*>(data)->specifications->front()),&(static_cast<tree_modelspec_struct*>(data)->values->at(0)));
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

double opt_rate_for_clades_function(const std::vector<double> &x, std::vector<double> &grad, void* data) {
    change_non_fixed(x,static_cast<tree_modelspec_struct*>(data)->values,static_cast<tree_modelspec_struct*>(data)->fixed);
    double LogLH = static_cast<tree_modelspec_struct*>(data)->tree->calculate_likelihood_rate_change_at_nodes(&(static_cast<tree_modelspec_struct*>(data)->specifications->front()),&(static_cast<tree_modelspec_struct*>(data)->values->at(0)), &(static_cast<tree_modelspec_struct*>(data)->values->at(static_cast<tree_modelspec_struct*>(data)->values->size()-static_cast<tree_modelspec_struct*>(data)->tree->n_taxon_sets())));
    return LogLH;

}


void change_non_fixed(const std::vector<double> &x, std::vector<double>* values, const std::set<unsigned int>* fixed) {
    /*#if DEBUG
        std::cerr << "No. fixed "<< fixed->size();
        std::cerr << " (";
        for (set<int>::iterator i=fixed->begin(); i!=fixed->end(); ++i) std::cerr << *i << ' ';
        std::cerr << ')' << endl;
        std::cerr << "No. parameters " << values->size() << endl;
        std::cerr << "Values to assign: ";
        for (vector<double>::const_iterator i=x.begin(); i!=x.end(); ++i) std::cerr << *i << ' ';
        std::cerr << endl;
    #endif //DEBUG*/
    unsigned int z(0);
    unsigned int j(0);
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

