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
#include <map>
#include "tree.h"
#include "decisiveness.h"
#include "file_parser.h"
//#include <vector>

using namespace std;

void help ();

struct tree_array {
    tree phylo;
    tree_array* next;
};

void pars_decisiveness_input( istream* input, string& genestring );

int main (int argc, char *argv []) {
    #if DEBUG
    cerr << "Debugging version!!!" << endl;
    #endif //DEBUG
    char method = 'c';
    float cut_off = 0.5;
    bool html = false;
    string database_name;
    ifstream database_file;
    string file_name;
    ifstream input_file;
    file_parser input(&cin);
    file_parser db_input(&cin);
    char print_format='w';
    //istream* input_stream = &std::cin;
    map <string,string> taxa_trans;
    string file_format;
    map <string,string> db_taxa_trans;
    string db_file_format;
    bool divide_by_n_trees(false);
    bool rooted(false);
    /////// Variables from superstat
    string genestring;
    int n_iterations(100);
    bool quiet(true);
    if (argc > 1) {
	/// Read arguments ///
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--compare")) {
                method = 'c';
                if ( i < argc-1 && argv[i+1][0] != '-') cut_off = atof(argv[++i]);
            }
            else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--robinson_foulds")) method = 'r';
            else if (!strcmp(argv[i],"-R") || !strcmp(argv[i],"--rooted")) rooted = true;
            else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--non_shared_tips")) method = 't';
	    else if (!strcmp(argv[i],"-a") || !strcmp(argv[i],"--add_to_support")) method = 'a';
	    else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--calculate_support")) { method = 'a'; divide_by_n_trees=true; }
	    //else if (!strcmp(argv[i],"-w") || !strcmp(argv[i],"--newick")) print_format = 'w';
	    //else if (!strcmp(argv[i],"-x") || !strcmp(argv[i],"--nexus")) print_format = 'x';
            else if (!strcmp(argv[i],"--html")) html = true;
            else if (!strcmp(argv[i],"--format")) {
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    string first; string second;
		    string* taking_it = &first;
		    for (unsigned int j=0; argv[i][j] != '\0'; ++j) {
			if (argv[i][j] == ',') taking_it = &second;
			*taking_it += argv[i][j];
		    }
		    if (!first.compare("nexus")) file_format = "nexus";
		    else if (!first.compare("newick")) file_format = "newick";
		    else {
			std::cerr << "Do not recognize format " << argv[i] << "." << endl;
			return 1;
		    }
		    if (second.empty()) db_file_format = file_format;
		    else if (!second.compare("nexus")) db_file_format = "nexus";
		    else if (!second.compare("newick")) db_file_format = "newick";
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
		    if (!strcmp(argv[i],"nexus") || !strcmp(argv[i],"nex") || (argv[i][0] == 'x' && argv[i][1] == '\0')) print_format = 'x';
		    else if (!strcmp(argv[i],"newick") || !strcmp(argv[i],"new") || (argv[i][0] == 'w' && argv[i][1] == '\0')) print_format = 'w';
		    else {
			std::cerr << "Do not recognize format " << argv[i] << "." << endl;
			return 1;
		    }
		}
		else std::cerr << "--output require nexus(nex or x) or newick (new or w) as additional argument" << endl;
	    }
	    else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--database")) {
		if ( i < argc-1 && (argv[i+1][0] != '-' || (argv[i+1][0] && argv[i+1][1] == '\0')))
		    database_name = argv[++i];
                else {
                    std::cerr << "-d/--database require a file name as next argument" << endl;
                    return 1;
                }
	    }
///// From superstat
	    else if (!strcmp(argv[i],"-D") || !strcmp(argv[i],"--decisiveness")) {
                method = 'D';
                if ( i < argc-1 && argv[i+1][0] != '-') genestring = argv[++i];
            }
	    else if (!strcmp(argv[i],"-i") || !strcmp(argv[i],"--iterations")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) n_iterations = atoi(argv[++i]);
		else {
		    cerr << "-i / --iterations must be followed by an integer number (e.g. -i 1000)." << endl;
		    return 1;
		}
	    }
    	    else if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose")) quiet = false;
////////////////
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
    if (!quiet) {
       std::cerr << "The program was called with the following command:" << endl;
       for (int i=0; i<argc; ++i) std::cerr << argv[i] << ' ';
       std::cerr << endl << endl;
    }
    /// Open files ///
    if (!file_name.empty()) {
        input_file.open(file_name.c_str(),std::ifstream::in);
        if (input_file.good())
            input.file_stream = &input_file;
        else {
            cerr << "Could not open file: " << file_name << endl;
            return 1;
        }
	if (!quiet) {
	    cerr << "Trees will be read from " << file_name << "." << endl;
	}
    }
    else if (!quiet) cerr << "Trees will be read from standard in." << endl;
    if (!file_format.empty()) {
	input.set_file_type(file_format.c_str());
    }
    else if (!input.set_file_type()) {
	if (!input.set_file_type("newick")) cerr << "Failed to set default file type." << endl;
    }
    if (!quiet) cerr << "Format of trees to compare is " << input.get_file_type() << "." << endl;
    if (!database_name.empty()) {
	if (!database_name.compare("-")) {
	    if (input.file_stream == &cin) {
		cerr << "Tree input must be from file if database is taken from standard in." << endl;
		return 1;
	    }
	    db_input.file_stream = &cin;
	    if (!quiet) {
		cerr << "Trees that will be compared against (database), will be read from stdin." << endl;
	    }
	}
	else {
	    database_file.open(database_name.c_str(),std::ifstream::in);
	    if (!database_file.good()) {
		cerr << "Could not open file: " << database_name << endl;
		return 1;
	    }
	    else db_input.file_stream = &database_file;
	    if (!quiet) {
		cerr << "Trees that will be compared against (database), will be read from " << database_name << "." << endl;
	    }
	}
	if (!db_file_format.empty()) {
    	    db_input.set_file_type(db_file_format.c_str());
	}
	else if (!db_input.set_file_type()) {
	    if (!db_input.set_file_type("newick")) cerr << "Failed to set default file type." << endl;
	}
	if (!quiet) cerr << "Format of trees to compare against (database) is " << db_input.get_file_type() << "." << endl;
    }
    /// Do calculation of decisiveness if that is what to be done
    if (method == 'D') {
	if (genestring.empty() && (input_file.good() || input.file_stream == &std::cin)) {
	    if (!quiet) cerr << "Parsing data for decisiveness estimation." << endl;
	    pars_decisiveness_input( input.file_stream, genestring );
	}
        if (genestring.empty()) {
            std::cerr << "A string of genes is needed to calculate the decisiveness (--decisiveness/-d), e.g. -d ITS,RPB2|ITS|ITS,RPB2." << endl;
            return 1;
        }
	#ifdef DEBUG
	cerr << genestring << endl;
	#endif //DEBUG
	if (!quiet) cerr << "Starting to estimate decisiveness." << endl; 
        decisiveness stat(&genestring);
        stat.on_random_tree( n_iterations );
        cout << "The gene sampling is decisive for " << stat.get_decisiveness() << " of the trees and " << stat.get_distinguished() << " of the branches." << endl;
	return 0;
    }
    /// Prepare to read nexus file
    char nexus_command = nexus_command::NON;
    if (input.test_file_type("nexus")) {
	if (input.move_to_next_X_block( nexus_block::TREES )) {
	    nexus_command = input.read_next_nexus_command();
	    if (nexus_command==nexus_command::TRANSLATE) {
		input.read_translate_parameters(taxa_trans);
		nexus_command = input.read_next_nexus_command();
	    }
	}
    }
    char db_nexus_command = nexus_command::NON;
    if (db_input.test_file_type("nexus")) {
	if (db_input.move_to_next_X_block( nexus_block::TREES )) {
	    db_nexus_command = db_input.read_next_nexus_command();
	    if (db_nexus_command==nexus_command::TRANSLATE) {
		db_input.read_translate_parameters(db_taxa_trans);
		db_nexus_command = db_input.read_next_nexus_command();
	    }
	}
    }
    /// Read trees
    tree_array* array = new tree_array;
    array->next = 0;
    tree_array* array_start = array;
    unsigned int read_trees(0);
    unsigned int print_trees(0);
    if (!quiet) cerr << "Reading trees ";
    while (1) {
	if (input.test_file_type("nexus") || input.test_file_type("newick")) {
    	    if (input.test_file_type("nexus")) {
		if (read_trees != 0) nexus_command = input.read_next_nexus_command();
		if (!(nexus_command==nexus_command::TREE && input.move_to_start_of_tree()))
		    break;
	    }
	    array->phylo.tree_file_parser( *(input.file_stream), taxa_trans, false );
	    ++read_trees;
	    if (array->phylo.empty()) break;
	    else {
		array->next = new tree_array;
		array = array->next;
		array->next = 0;
	    }
	}
	else { cerr << "Unrecognized file format." << endl; return 1; }
	if (!quiet) cerr << '.';
    }
    if (!quiet) cerr << endl << "Starting to compare trees." << endl;
    array = array_start;
    /// Beginning messages for the different methods
    if (method == 'c' && html) {
	std::cout << "<!DOCTYPE html>" <<  endl << "<html>" << endl << "<body>" << endl;
	std::cout << "<h1>Conflicts between trees</h1>" << endl;
	std::cout << "<p>Tips in green are supported as monophyletic in the first tree, while tips in red are supported as nested within the green in the second tree.</p>" << endl;
    }
    /// Loop over trees
    int i=0;
    unsigned int read_db_trees(0);
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
	if (db_input.test_file_type("nexus")) {
	    if (db_input.move_to_next_X_block( nexus_block::TREES )) {
		db_nexus_command = db_input.read_next_nexus_command();
		if (db_nexus_command==nexus_command::TRANSLATE) {
		    //db_input.read_translate_parameters(db_taxa_trans);
		    db_nexus_command = db_input.read_next_nexus_command();
		}
	    }
	}
	/// Loop over trees to compare to
	while (1) {
	    ++j;
	    convert.clear(); // clear to-string-converter for new input
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
		if (db_input.test_file_type("nexus")) {
		    #ifdef DEBUG
		    std::cerr << "It is in nexus format. " << db_taxa_trans.size() << ' ' << read_db_trees << std::endl;
		    #endif //DEBUG
		    if (read_db_trees != 0) db_nexus_command = db_input.read_next_nexus_command();
		    if (!(db_nexus_command==nexus_command::TREE && db_input.move_to_start_of_tree()))
			break;
		    #ifdef DEBUG
		    cerr << db_nexus_command << endl;
		    #endif //DEBUG
		}
		tree_j = new tree;
		tree_j->tree_file_parser( *(db_input.file_stream), db_taxa_trans, false );
		++read_db_trees;
		#ifdef DEBUG
		std::cerr << "Read tree with: " << tree_j->n_tips() << std::endl;
		#endif //DEBUG
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
		std::cout << "Robinson-Foulds metric between tree " << i << " and tree " << j << " (per internal node): ";
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
	    else if (method == 'a') array->phylo.add_to_support(tree_j, rooted);
	    if (database_name.empty()) comp_tree = comp_tree->next;
	    else delete tree_j;
	}
	if (method == 'r') std::cout << "Sum of Robinson-Foulds metric for tree " << i << " (per internal node, number of comparisons): " << sum << " (" << normalized_sum << ", " << n_comp << ')' << std::endl;
	if (method == 'a') {
	    if (database_name.empty()) {
		array->phylo.add_one_to_each_support(rooted);
		//--j;
	    }
	    #ifdef DEBUG
	    cerr << "N trees compared: " << j << endl;
	    #endif //DEBUG
	    if (divide_by_n_trees) array->phylo.divide_support_by(float(j-1));
	    if (print_format == 'w') array->phylo.print_newick();
	    else if (print_format == 'x') {
		if (print_trees == 0) array->phylo.print_nexus_tree_intro(taxa_trans);
		stringstream ss;
		ss << "tree" << read_trees;
		ss << '_' << read_trees;
		string tree_name(ss.str());
		array->phylo.print_tree_to_nexus( tree_name, true, true, taxa_trans );
	    }
	    ++print_trees;
	}
	array = array->next;
    }
    if (method == 'a' && print_format == 'x') cout << "End;" << endl;
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
    std::cout << "Contree " << VERSION << " is a command line program for comparing trees. It can compare" << endl;
    std::cout << "trees in one file/input against each other or compare the trees in one file/"<< endl;
    std::cout << "input to the trees in another file/input." << endl;
    std::cout << "(c) Martin Ryberg " << YEAR << "." << endl << endl;
    std::cout << "Usage:" << endl << "contree [arguments] file.trees" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--add_to_support / -a           add one to the value of the internal node for" << endl;
    std::cout << "                                each tree that that split is present in." <<endl;
    std::cout << "--calculate_support / -s        same as --add_to_support but divide the number" << endl;
    std::cout << "                                by number of trees compared against." << endl;
    std::cout << "--compare / -c [value]          output conflicting splits where at least one" << endl;
    std::cout << "                                branch support the conflict with more than given" << endl;
    std::cout << "                                support, e.g. -c 0.7." << endl;
    std::cout << "--database / -d [file]          give a second file of trees to compare against" << endl;
    std::cout << "                                instead of comparing within the ordinary input." << endl;
    std::cout << "                                If - is given input will be taken from standard" << endl;
    std::cout << "                                in instead of file, e.g. -d -. In that case the" << endl;
    std::cout << "                                input compared to must  be given in a file (by" << endl;
    std::cout << "                                -f or as last argument)." << endl;
    std::cout << "--decisiveness / -D [genes]     calculates proportion of random trees for which" << endl;
    std::cout << "                                given gene sampling is decisive and the mean" << endl;
    std::cout << "                                proportion of branches that are distinguished." << endl;
    std::cout << "                                (Sandersson et al. 2010, BMC Evol. Biol. 19:" << endl;
    std::cout << "                                155). The genes for each taxon are given as a" << endl;
    std::cout << "                                comma (,) separated string, the genes for" << endl;
    std::cout << "                                different taxon are separated by a bar (|). E.g" << endl;
    std::cout << "                                -D 'ITS,RPB2|ITS|ITS,RPB2|RPB2|RPB2|ITS', or in" << endl;
    std::cout << "                                a file (given by -f or as last argument) with a" << endl;
    std::cout << "                                comma separated string with the genes for each" << endl;
    std::cout << "                                taxa on a separate row." << endl;
    std::cout << "--iterations / -i [number]      give numbers of iterations to do when" << endl;
    std::cout << "                                calculating decisiveness, e.g. -i 1000." << endl; 
    std::cout << "--file / -f [file]              give file name for trees, or decisiveness, e.g." << endl;
    std::cout << "                                -f file.tree." << endl;
    std::cout << "--format [newick/nexus]         give format of input, e.g. --format nexus. If no" << endl;
    std::cout << "                                format is given and the input is a file contree" << endl;
    std::cout << "                                will try to guess the format, if it is through" << endl;
    std::cout << "                                standard in it will assume newick format. A" << endl;
    std::cout << "                                separate file format can be given for the" << endl;
    std::cout << "                                database file after a ',', e.g. --format" << endl;
    std::cout << "                                newick,nexus. If only one format is given, it" << endl;
    std::cout << "                                will be used for both tree files/inputs." << endl;
    std::cout << "--help / -h                     print this help." << endl;
    std::cout << "--html                          give output as tree in html (svg) format with" << endl;
    std::cout << "                                conflicting tips coloured green and red when" << endl;
    std::cout << "                                doing --compare." << endl;
    std::cout << "--non_shared_tips / -t          print tip names not present in the other tree." << endl;
    std::cout << "--output [newick/nexus]         give tree format for output, nexus (nex or x for" << endl;
    std::cout << "                                short) or newick (new or w for short), e.g" << endl;
    std::cout << "                                --output x. (default w)." << endl;
    std::cout << "--robinson_foulds / -r          compute Robinson-Foulds metric between trees." << endl;
    std::cout << "--rooted / -R                   treat trees as rooted for --add_to_support and" << endl; 
    std::cout << "                                --calculate_support." << endl;
    std::cout << "--verbose / -v                  get additional output." << endl;
    std::cout << endl;
}

void pars_decisiveness_input( istream* input, string& genestring ) {
    while (*input) {
	string line;
	*input >> line;
	if (!line.empty()) {
	    if (!genestring.empty()) genestring += '|';
	    genestring += line;
	}
    }
}
