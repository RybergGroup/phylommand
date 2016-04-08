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
#include "tree.h"
#include "nj_tree.h"
#include "character_vector.h"
#include "matrix_parser.h"

using namespace std;

void help ();

int main (int argc, char *argv []) {
    bool lables = true;
    char method = 'p';
    bool random = false;
    bool print_br_length(true);
    string tree_file_name;
    ifstream tree_file;
    istream* tree_stream = &std::cin;
    string data_file_name;
    ifstream data_file;
    istream* data_stream = &std::cin;
    string alphabet_file_name;
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
            else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--likelihood")) {
		method = 'l';
	    }
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
    if (method == 'p' || method == 's') {
	vector<character_vector> characters;
////////////////////
	map<char, bitset<SIZE> > alphabet;
	set_alphabet_binary(alphabet);
	matrix_parser data_parser(*data_stream, characters, alphabet);
	data_parser.pars();
	#ifdef DEBUG
	std::cerr << "Max number of characters: " << characters.begin()->max_n_char() << std::endl;
	#endif //DEBUG
        if (data_file.is_open()) data_file.close();
/////////////////////////////////
	if (method == 'p') {
	    while (1) {
		tree tree;
		tree.tree_file_parser( *tree_stream );
		if (tree.empty()) break;
		std::cout << tree.fitch_parsimony( characters, print_br_length ) << std::endl;
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
    }
    if (method == 'l') {
	std::cout << "Sorry, likelihood is not available yet." << std::endl;
	return 0;
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
    std::cout << "--data_file / -d [file name]          give the data file"<< std::endl;
    std::cout << "--neighbour_joining / -n              compute neighbour joining tree for given data. The data should be" << std::endl;
    std::cout << "                                           a left triangular similarity matrix." << std::endl;
    std::cout << "--no_branch_length / -0               do not print branch lengths and do not calculate branch lengths for" << std::endl;
    std::cout << "                                           parsimony trees" << endl;
    std::cout << "--parsimony / -p                      calculate parsimony score for given tree and data." << std::endl;
    std::cout << "--likelihood / -l (not yet valid)" << std::endl;
    std::cout << "--step_wise / -s                      do parsimony stepwise addition." << std::endl;
    std::cout << "--tree_file / -t [file name]          give tree file name" << std::endl;
    std::cout << "--alphabet_file / -a (not yet valid)" << std::endl;
    std::cout << "--help / -h                           print this help." << endl;
    std::cout << "--no_lable / -l                       will tell treeator that there are no taxon labels in the matrix." << endl;
    std::cout << "--random / -r                         do stepwise addition in random order." << endl;
    std::cout << endl;
}

