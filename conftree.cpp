/********************************************************************
Copyright (C) 2014 Martin Ryberg

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
#include "tree.h"
//#include <vector>

using namespace std;

void help ();

struct tree_array {
    tree phylo;
    tree_array* next;
};

int main (int argc, char *argv []) {
    char method = 'c';
    float cut_off = 0.5;
    bool html = false;
    string database_name;
    ifstream database_file;
    string file_name;
    ifstream input_file;
    istream* input_stream = &std::cin;
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--compare")) {
                method = 'c';
                if ( i < argc-1 && argv[i+1][0] != '-') cut_off = atof(argv[++i]);
            }
            else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--Robinson–Foulds")) method = 'r';
            else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--non_shared_tips")) method = 't';
	    else if (!strcmp(argv[i],"-a") || !strcmp(argv[i],"--add_to_support")) method = 'a';
            else if (!strcmp(argv[i],"--html")) html = true;
	    else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--database")) {
		if ( i < argc-1 && argv[i+1][0] != '-' )
		    database_name = argv[++i];
                else {
                    std::cerr << "-d/--database require a file name as next argument" << endl;
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
    if (!database_name.empty()) {
        database_file.open(database_name.c_str(),std::ifstream::in);
        if (!database_file.good()) {
            cerr << "Could not open file: " << database_name << endl;
            return 1;
        }
    }
    // Read trees
    tree_array* array = new tree_array;
    array->next = 0;
    tree_array* array_start = array;
    while (1) {
        array->phylo.tree_file_parser( *input_stream );
        if (array->phylo.n_tips() < 2) break;
        else {
            array->next = new tree_array;
            array = array->next;
            array->next = 0;
        }
    }
    array = array_start;
    // Beginning messages for the different methods
    if (method == 'c' && html) {
	std::cout << "<!DOCTYPE html>" <<  endl << "<html>" << endl << "<body>" << endl;
	std::cout << "<h1>Conflicts between trees</h1>" << endl;
	std::cout << "<p>Tips in green are supported as monophyletic in the first tree, while tips in red are supported as nested within the green in the second tree.</p>" << endl;
    }
    //Loop over trees
    int i=0;
    while (array != 0 && array->phylo.n_tips() > 1) {
	unsigned int sum = 0;
	double normalized_sum = 0.0;
	++i;
	std::stringstream convert;
	convert << i;
	string name_i = convert.str();
	tree_array* comp_tree = array_start;
	int j=0;
	int n_comp = 0;
	database_file.clear();
	database_file.seekg(0,database_file.beg);
	// Loop over trees to compare to
	while (1) {
	    ++j;
	    convert.clear();
	    convert.str(std::string());
	    convert << j;
	    string name_j = convert.str();
	    tree* tree_j;
	    if (database_name.empty()) {
		if (comp_tree == 0 || comp_tree->phylo.empty()) break;
		if (comp_tree == array) {
		    comp_tree = comp_tree->next;
		    continue;
		}
		tree_j = &(comp_tree->phylo);
	    }
	    else {
		#ifdef DEBUG
		std::cerr << "Database file given." << std::endl;
		#endif //DEBUG
		tree_j = new tree;
		tree_j->tree_file_parser( database_file );
		if (tree_j->empty()) {
		    delete tree_j;
		    break;
		}
	    }
	    if (method == 'c') {
		if (html) std::cout << "<h2>" << endl;
		std::cout << "Checking conflicts between tree " << i << " and tree " << j << ":" << endl;
		if (html) {
		    std::cout << "</h2>" << endl;
		}
		array->phylo.print_conflict_clades_reduced( tree_j, cut_off, html );
	    }
	    else if (method == 'r') {
		std::cout << "Robinson–Foulds metric between tree " << i << " and tree " << j << " (per internal node): ";
		unsigned int n_shared_taxa = array->phylo.shared_tips( tree_j );
		if (n_shared_taxa == 0) std::cout << "0 (No shared taxa)" << endl;
		if (n_shared_taxa < 4) std::cout << "0 (Less than 4 taxa in common)" << endl;
		else {
		    unsigned int rf = array->phylo.robinson_fould( tree_j );
		    double normalized = rf/((static_cast<double>(n_shared_taxa)-3)*2);
		    std::cout << rf << " (" << normalized << ')' << std::endl;
		    sum += rf;
		    normalized_sum += normalized;
		    ++n_comp;
		}
	    }
	    else if (method == 't')
		array->phylo.tips_not_shared( tree_j, name_i, name_j );
	    else if (method == 'a') array->phylo.add_to_support(tree_j);
	    if (database_name.empty()) comp_tree = comp_tree->next;
	    else delete tree_j;
	}
	if (method == 'r') std::cout << "Sum of Robinson–Foulds metric for tree " << i << " (per internal node, and number of comparisons): " << sum << " (" << normalized_sum << ", " << n_comp << ')' << std::endl;
	if (method == 'a') array->phylo.print_newick();
	array = array->next;
    }
    if (html) std::cout << "</body>" << endl << "</html>" << endl;
    array=array_start;
    tree_array* array_next;
    while (array->next != 0) {
        array_next = array->next;
        delete array;
        array = array_next;
    }
}

void help () {
    std::cout << "Conftree is a command line program for comparing trees." << endl;
    std::cout << "The program take two trees in newick format as indata through standard in." << endl;
    std::cout << "(c) Martin Ryberg 2015." << endl << endl;
    std::cout << "Usage:" << endl << "conftree [arguments] < file.trees" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--add_to_support / -a                                      add one to the value of internal node for each tree the split is present in." <<endl;
    std::cout << "--compare / -c [float value]                               output conflicting splits where at least one branch support the conflict with more than given support, e.g. -f 0.7." <<endl;
    std::cout << "--database / -d [file name]                                give a second file of trees to compare agains instead of comparing within the ordinary input." << endl;
    std::cout << "--file/-f [file name]                                      give file name for trees, e.g. -f file.tree." << endl;
    std::cout << "--help / -h                                                print this help." << endl;
    std::cout << "--html                                                     give output as tree in html (svg) format with conflicting tips coloured green and red." << endl;
    std::cout << "--non_shared_tips / -t                                     print tip names not present in other tree." << endl;
    std::cout << "--Robinson–Foulds / -r                                     compute Robinson–Foulds metric between trees." << endl;
    std::cout << endl;
}

