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
#include <map>
#include <vector>
#include <limits.h>
#include "tree.h"
#include "file_parser.h"
//#include "supmar/support_functions.h"

using namespace std;

void help ();
string pars_ARGV_string ( char input [] );

int main (int argc, char *argv []) {
    #ifdef DEBUG
    cerr << "Debug mode" << endl;
    #endif //DEBUG
    char method('0');
    char tree_source('s');
    char output_format('w');
    char inverse_taxa('n');
    unsigned int tree_interval_start = 0;
    unsigned int tree_interval_end = UINT_MAX;
    string taxastring;
    vector<string> taxon_vector;
    string separator(",");
    string separator2(" ");
    bool print_scalebar(true);
    bool print_nodelabel(true);
    // svg
    float svg_width(-1.0);
    float svg_height(-1.0);
    float textoffset(5.0);
    int strokewidth(1);
    int fontsize(5);
    string font("Arial");
    string tip_colors;
    /////
    string stats_param;
    bool read_figtree_annotations = false;
    bool print_br_length(true);
    bool flag(false);
    bool right_ladderize(true);
    bool read_phylomand_block = false;
    float value(1.0);
    float cut_off(0.0);
    int n_taxa(0);
    unsigned int n_rand_trees(1);
    unsigned int branch_no(0);
    map<string, vector<char> > matrix;
    string file_name;
    //char input_format('w');
    string file_format;
    ifstream input_file;
    //istream* input_stream = &std::cin;
    file_parser input(&cin);
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--mid_point_root")) method = 'm';
            else if (!strcmp(argv[i],"-o") || !strcmp(argv[i],"--outgroup_root")) {
                method = 'g';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = pars_ARGV_string( argv[++i]);
		else { cerr << "-o/--outgroup_root need an outgroup as next argument" << endl; return 1; }
            }
	    else if (!strcmp(argv[i],"--relaxed_outgroup_root")) {
                method = 'G';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = pars_ARGV_string( argv[++i]);
                else { cerr << "--relaxed_outgroup_root need an outgroup as next argument" << endl; return 1; }
            }
            else if (!strcmp(argv[i],"--get_relaxed_outgroup")) {
                method = 'H';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = pars_ARGV_string( argv[++i]);
		else { cerr << "--get_relaxed_outgroup need an outgroup as next argument" << endl; return 1; }
            }
            else if (!strcmp(argv[i],"--is_monophyletic")) {
                method = 'M';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = pars_ARGV_string( argv[++i]); //argv[++i];
            }
            else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--ladderize")) {
                method = 'l';
                if ( i < argc-1 && argv[i+1][0] != '-') {
                    if (!(strcmp(argv[i+1],"l")) or !(strcmp(argv[i+1],"L")) or !(strcmp(argv[i+1],"left"))) { right_ladderize = 0; ++i; }
                    else if (!(strcmp(argv[i+1],"r")) or !(strcmp(argv[i+1],"R")) or !(strcmp(argv[i+1],"right"))) ++i;
                }
            }
            else if (!strcmp(argv[i],"-a") || !strcmp(argv[i],"--branch_lengths")) {
                method = 'a';
                if ( i < argc-1 && argv[i+1][0] != '-') {
                    separator = argv[++i];
                    if (i < argc-1 && argv[i+1][0] != '-') {
                        if (argv[i+1][0] == 'r' && strlen(argv[i+1])==1) { flag=true; ++i; }
                        else if (argv[i+1][0] == 'n' && strlen(argv[i+1])==1) { flag=false; ++i; }
                    }
                }
            }
            else if (!strcmp(argv[i],"-z") || !strcmp(argv[i],"--distances_to_root")) {
                method = 'o';
                if ( i < argc-1 && argv[i+1][0] != '-') {
                    separator = argv[++i];
                    if ( i < argc-1 && argv[i+1][0] != '-') {
                        separator2 = argv[++i];
                    }
                    else separator2 = "\t";
                }
                else separator = "\n";
            }
            else if (!strcmp(argv[i],"-p") || !strcmp(argv[i],"--patristic_distances")) {
                method = 'p';
                if ( i < argc-1 && argv[i+1][0] != '-') {
                    separator = argv[++i];
                    if ( i < argc-1 && argv[i+1][0] != '-') {
                        separator2 = argv[++i];
                    }
                    else separator2 = "\t";
                }
                else separator = "\n";
            }
            else if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"--set_branch_lengths")) {
                method = 'b';
                if ( i < argc-1 && argv[i+1][0] != '-') value = atof(argv[++i]);
            }
            else if (!strcmp(argv[i],"-u") || !strcmp(argv[i],"--multiply_branch_lengths")) {
                method = 'u';
                if ( i < argc-1 && argv[i+1][0] != '-') value = atof(argv[++i]);
            }
	    else if (!strcmp(argv[i],"-U") || !strcmp(argv[i],"--multiply_branch_lengths_until")) {
		method = 'U';
		if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    string temp;
		    for (unsigned int j=0; argv[i][j] != '\0'; ++j) {
			if (argv[i][j] == ',') {
			    cut_off = atof(temp.c_str());
			    temp.clear();
			}
			else temp += argv[i][j];
		    }
		    value = atof(temp.c_str());
		}
		//cerr << cut_off << ' ' << value << endl;
	    }
	    else if (!strcmp(argv[i],"-V") || !strcmp(argv[i],"--multiply_branch_lengths_clade")) {
                method = 'V';
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    string temp;
		    char mode='s';
		    for (unsigned int j=0; argv[i][j] != '\0'; ++j) {
			if (mode == 's' && (argv[i][j] == ':' || argv[i][j] == ';')) { 
			    value = atof(temp.c_str());
			    mode = 't';
			    temp.clear();
			}
			else if (mode == 't' && (argv[i][j] == ':' || argv[i][j] == ';')) {
			    taxon_vector.push_back(temp);
			    temp.clear();
			}
			else temp += argv[i][j];
		    }
		    if (mode == 't' && !temp.empty()) taxon_vector.push_back(temp);
		    else if (mode == 's' && !temp.empty()) value = atof(temp.c_str());
		    else {
			cerr << "Parsing error reading argument to --multiply_branch_lengths_clade / -V." << endl;
			return 1;
		    }
		}
	    }
            else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--get_tip_names")) {
                method = 't';
                if ( i < argc-1 && argv[i+1][0] != '-') separator = argv[++i];
            }
            else if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--change_names")) {
                method = 'c';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = pars_ARGV_string( argv[++i]); //argv[++i];
            }
            else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--drop_tips")) {
                method = 'd';
                if ( i < argc-1 && argv[i+1][0] != '-' ) taxastring = pars_ARGV_string( argv[++i]); //argv[++i];
            }
            else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--random_tree")) {
                tree_source = 'r';
                if ( i < argc-1 && argv[i+1][0] != '-' ) n_taxa = atoi(argv[++i]);
            }
            else if (!strcmp(argv[i],"--interval")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    string argument(argv[++i]);
		    size_t separator = argument.find("-");
		    tree_interval_start = atoi(argument.substr(0,separator).c_str());
		    if (separator < argument.length()) tree_interval_end = atoi(argument.substr(separator+1).c_str());
		    if (tree_interval_end < tree_interval_start) {
			cerr << "The end of the interval (--interval) has to be smaler than the begining." << endl;
			return 1;
		    }
		}
		else {
		    cerr << "--interval require at least one number to set the interval to print." << endl;
		    return 1;
		}
            }
            else if (!strcmp(argv[i],"--n_supported")) {
                method='Y';
                if ( i < argc-1 && argv[i+1][0] != '-') value = atof(argv[++i]);
            }
            else if (!strcmp(argv[i],"--aMPL")) {
                method = 'C';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = pars_ARGV_string( argv[++i]); //argv[++i];
            }
            else if (!strcmp(argv[i],"--nni")) {
                method = 'I';
		if (i+1 < argc && argv[i+1][0] != '-') {
		    ++i;
		    if (!strcmp(argv[i],"all")) flag = true;
		    else branch_no = atoi(argv[i]);
		}
	    }
            else if (!strcmp(argv[i],"--clocklikness")) method = 'K';
	    else if (!strcmp(argv[i],"--read_figtree_annotations")) read_figtree_annotations = true;
            else if (!strcmp(argv[i],"--matrix_representation")) method = 'R';
            else if (!strcmp(argv[i],"--internal_node_stats")) {
		method = 'T';
		if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    int j=0;
		    while (argv[i][j] != '\0' && argv[i][j] != '/') ++j;
		    string temp;
		    for (int k=0; k<j; ++k) temp+=argv[i][k];
		    value = atof(temp.c_str());
		    if (argv[i][j] != '\0') ++j;
		    while (argv[i][j] != '\0') stats_param += argv[i][j++];
		}
	    }
            else if (!strcmp(argv[i],"-n") || !strcmp(argv[i],"--number_of_taxa")) method = 'n';
            else if (!strcmp(argv[i],"-D") || !strcmp(argv[i],"--depth")) method = 'D';
            else if (!strcmp(argv[i],"-0") || !strcmp(argv[i],"--no_branch_length")) print_br_length = false;
            else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--sum_branch_length")) method = 's';
            else if (!strcmp(argv[i],"-i") || !strcmp(argv[i],"--inverse")) inverse_taxa = 'y';
            else if (!strcmp(argv[i],"-w") || !strcmp(argv[i],"--newick")) output_format = 'w';
            else if (!strcmp(argv[i],"-x") || !strcmp(argv[i],"--nexus")) output_format = 'x';
            else if (!strcmp(argv[i],"--get_branch_numbers")) method = 'N';
            else if (!strcmp(argv[i],"--format")) {
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    if (!strcmp(argv[i],"nexus")) file_format = "nexus";
		    else if (!strcmp(argv[i],"newick")) file_format = "newick";
		    else {
			std::cerr << "Do not recognize format " << argv[i] << "." << endl;
			return 1;
		    }
		}
		else std::cerr << "--format require nexus or newick as additional argument" << endl;
	    }
            else if (!strcmp(argv[i],"-g") || !strcmp(argv[i],"--svg")) {
		output_format = 's';
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    unsigned int j=0;
		    string temp;
		    string param;
		    while (1) {
		    	if (argv[i][j] == '&' || argv[i][j] == '\0') {
			    if (!param.compare("width")) svg_width = atof(temp.c_str());
			    else if (!param.compare("height")) svg_height = atof(temp.c_str());
			    else if (!param.compare("offset")) textoffset = atof(temp.c_str());
			    else if (!param.compare("strokewidth"))  strokewidth = atoi(temp.c_str());
			    else if (!param.compare("fontsize")) fontsize = atoi(temp.c_str());
			    else if (!param.compare("font")) font = temp;
			    else if (!param.compare("tipcolor") || !param.compare("tipcolour")) tip_colors = temp;
			    else { cerr << "Unrecognized parameter for svg output: " << param << endl; return 1; }
			    if (argv[i][j] == '\0') break;
			    temp.clear();
			    param.clear();
			}
			else if (argv[i][j] == ':') {
			    param = temp;
			    temp.clear();
			}
			else temp += argv[i][j];
			++j;
		    }
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
    if ( (method == 'c' || method == 'd' || method == 'M') && taxastring.empty() ) {
        std::cerr << "Selected method reqire a taxon string. Use -h to print help." << endl;
        return 1;
    }
    #ifdef DEBUG
    if (!taxastring.empty()) {
	cout << "Taxon string: " << taxastring << endl;
    }
    #endif //DEBUG
    if ( !separator.compare("\\n") ) separator = "\n";
    else if ( !separator.compare("\\t") ) separator = "\t";
    else if ( !separator.compare("\\r") ) separator = "\r";
    if (!file_name.empty()) {
	#ifdef DEBUG
	cerr << "File name " << file_name << endl;
	#endif //DEBUG
	/*if (!input_file.open_file(file_name)) {
	    cerr << "Could not open file: " << file_name << endl;
	    return 1;
	}*/
	input_file.open(file_name.c_str(),std::ifstream::in);
	if (input_file.good())
	    input.file_stream = &input_file;
	else {
	    cerr << "Could not open file: " << file_name << endl;
	    return 1;
	}
    }
    bool print_tree = true;
    map <string,string> taxa_trans;
    //if (input_stream != &std::cin && file_format.empty())
	// file_format = support_functions::get_file_type( *input_stream );
    #ifdef DEBUG
    cerr << "File type " << input.get_file_type() << endl;
    #endif //DEBUG
    if (!file_format.empty())  {
	input.set_file_type(file_format.c_str());
	#ifdef DEBUG
	cerr << "Try to set file type " << file_format << endl;
	#endif //DEBUG
    }
    else if (!input.set_file_type()) {
	#ifdef DEBUG
	cerr << "Try to set default file type" << endl;
	#endif //DEBUG
	if (!input.set_file_type("newick")) cerr << "Failed to set default file type." << endl;
	
    }
    #ifdef DEBUG
    cerr << "File type " << input.get_file_type() << endl;
    #endif //DEBUG
    if (method=='C' && !taxastring.compare("nexus")) { // If doing aMPL and nexus is given
	read_phylomand_block = true;
	taxastring.clear();
    }
    if (read_phylomand_block) {
	#ifdef DEBUG
	cerr << "Trying to read phylommand block." << endl;
	#endif //DEBUG
	if (input.test_file_type("nexus")) {
    	    if (input.move_to_next_X_block( nexus_block::PHYLOMMAND )) {
		char command('S');
		map<string,string> taxon_set;
		map<string,string> ages;
		while (command != nexus_command::NON) {
		    command = input.read_next_nexus_command();
		    #ifdef DEBUG
		    cerr << command << endl;
		    #endif //DEBUG
		    if (command == nexus_command::MRCA) {
			if (!input.read_mrca(taxon_set)) cerr << "Failed to read a taxset command in PHYLOMMAND block." << endl;
		    }
		    else if (command == nexus_command::FIXAGE) {
			if (!input.read_fixage(ages)) cerr << "Failed to read a fixage in PHYLOMMAND block." << endl;
		    }
		    else if (command == nexus_command::END) {
			#ifdef DEBUG
			cerr << "Found end of phylommand block." << endl;
			#endif //DEBUG
			break;
		    }
		    else input.move_to_end_of_command();
		}
		if (method=='C') {
		    #ifdef DEBUG
	    	    cerr << "N fixed nodes = " << ages.size() << " N taxon sets = " << taxon_set.size() << endl;
    		    #endif //DEBUG
		    for (map<string,string>::const_iterator i=ages.begin(); i!=ages.end(); ++i) {
			map<string,string>::const_iterator taxa = taxon_set.find(i->first);
			if (!taxastring.empty()) taxastring += ';';
			if (taxa != taxon_set.end())
			    taxastring += taxa->second;
			else
			    taxastring += i->first;
			taxastring += ':';
			taxastring += i->second;
		    }
		}
	    }
	    else {
		cerr << "Did not find any phylommand block." << endl;
		return 1;
	    }
	}
	else {
	    cerr << "File not in nexus format, so can not read PHYLOMMAND block from file." << endl;
	    return 1;
	}
    }
    char nexus_command = nexus_command::NON;
    if (input.test_file_type("nexus")) { //!file_format.compare("nexus")) {
	if (input.move_to_next_X_block( nexus_block::TREES )) {
	    nexus_command = input.read_next_nexus_command();
	    if (nexus_command==nexus_command::TRANSLATE) {
		input.read_translate_parameters(taxa_trans);
		nexus_command = input.read_next_nexus_command();
	    }
	}
    }
    unsigned int read_trees = 0;
    while (1) {
	vector<tree> in_tree;
	in_tree.push_back(tree());
	// Get tree
	if (tree_source == 's') {
	    //if (input_format == 'w' || input_format == 'x') {
	    if (input.test_file_type("nexus") || input.test_file_type("newick")) {
		//if (input_format == 'x') {
		if (input.test_file_type("nexus")) {
		    if (read_trees != 0) nexus_command = input.read_next_nexus_command();//nexus_support::read_tree_block_command(*input_stream);
		    //if (!(nexus_command==nexus_support::TREE && nexus_support::move_to_start_of_tree(*input_stream)))
		    if (!(nexus_command==nexus_command::TREE && input.move_to_start_of_tree()))
			break;
		}
		in_tree.back().tree_file_parser( *(input.file_stream), taxa_trans, read_figtree_annotations );
		//in_tree.back().tree_file_parser( *input_stream, taxa_trans, read_figtree_annotations );
		++read_trees;
	    }
	    else { std::cerr << "Do not recognize tree format" << std::endl; return 1; }
	}
	else if (tree_source == 'r') {
            if (n_rand_trees<1) break;
	    in_tree.back().rand_topology(n_taxa);
	    --n_rand_trees;
	    ++read_trees;
	}
	if (in_tree.back().n_tips() < 2) break;
	// Make preparations
	if (inverse_taxa == 'y') taxastring = in_tree.back().not_present( taxastring );
    
	// Do manipulation
	if (method == 'm') in_tree.back().midpoint_root( );
	else if (method == 'M') {
	    if ( in_tree.back().is_monophyletic( taxastring ) ) std::cout << "Yes" << endl;
	    else std::cout << "No" << endl;
	    print_tree = false;
	}
	else if (method == 't') {
	    std::cout << in_tree.back().tip_names( separator ) << endl;
	    print_tree = false;
	}
	else if (method == 'g') in_tree.back().outgroup_root( taxastring );//std::cerr << "No function for outgroup rooting yet available" << endl; //tree.name_clust_cout( name_position, separator);
	else if (method == 'G') in_tree.back().relaxed_outgroup_root( taxastring );//std::cerr << "No function for outgroup rooting yet available" << endl; //tree.name_clust_cout( name_position, separator);
	else if (method == 'H') {
	    in_tree.back().print_relaxed_outgroup( taxastring );
	    cout << endl;
	    print_tree = false;
	}
	else if (method == 'o') {
	    in_tree.back().print_distance_to_root_for_all ( separator2, separator );
	    print_tree = false;
	}
	else if (method == 'l') in_tree.back().ladderize(right_ladderize);
	else if (method == 'd') in_tree.back().drop_tips( taxastring );
	else if (method == 'c') in_tree.back().change_tip_names( taxastring );
	else if (method == 'n') {
	    std::cout << in_tree.back().n_tips() << endl;
	    print_tree = false;
	}
	else if (method == 'D') {
	    std::cout << in_tree.back().longest_to_tip() << endl;
            print_tree = false;
	}
	else if (method == 's') {
	    std::cout << in_tree.back().sum_br_length() << endl;
	    print_tree = false;
	}
	else if (method == 'u') in_tree.back().multiply_br_length( value );
	else if (method == 'U') in_tree.back().multiply_br_length_cut_off( cut_off, value );
	else if (method == 'V') in_tree.back().multiply_br_length_clades( taxon_vector, value );
	else if (method == 'b') in_tree.back().set_br_length( value );
	else if (method == 'a') {
	    in_tree.back().print_branch_lengths ( separator, flag );
	    print_tree = false;
	}
	else if (method == 'I') {
	    unsigned int start = 1;
	    unsigned int n_tips = in_tree.back().n_tips();
	    if (n_tips < 4) std::cerr << "Warning!!! Can not do NNI on tree with fewer than 4 tips." << std::endl;
	    if (branch_no > (n_tips-3))  std::cerr << "Warning!!! Branch to swap on is out of bound (" << branch_no << " > " << (n_tips-3) << ")." << std::endl;
	    if (flag) branch_no = (n_tips-3);
	    else start = branch_no;
	    for (unsigned int i=start; i<=branch_no; ++i) {
		tree swapOne;
		tree swapTwo;
		in_tree[0].nni(i,swapOne,swapTwo);
		in_tree.push_back(swapOne);// swap_one;
		in_tree.push_back(swapTwo);// swap_one;
	    }
	}
	else if (method == 'p') {
	    in_tree.back().print_distances_to_all_for_each_taxa ( separator2, separator );
	    print_tree = false;
	}
	else if (method == 'Y') {
	    std::cout << in_tree.back().n_supported(value) << endl;
	    print_tree = false;
	}
       	else if (method == 'C') {
	    if (taxastring.empty()) taxastring="root:1";
	    #ifdef DEBUG
	    cerr << taxastring << endl;
	    #endif //DEBUG
	    in_tree.back().adjustedMPL(taxastring);
	}
       	else if (method == 'K') in_tree.back().test_clock_likness();
	else if (method == 'T') {
	    in_tree.back().internal_nodes_stat(stats_param,true,true,true,true,true,true,value,true,true,true,true);
	    print_tree = false;
	}
	else if (method == 'R') {
	    in_tree.back().add_to_matrix_representation(matrix);
	    print_tree = false;
	    unsigned int length_first = matrix.begin()->second.size();
	    bool warn = false;
	    for (map<string, vector<char> >::const_iterator i = matrix.begin(); i!= matrix.end(); ++i) {
		if (i->second.size() != length_first) warn = true;
	    }
	    if (warn) std::cerr << "WARNING!!! Matrix rows of different lengths!!!" << std::endl;
	}
	else if (method == 'N') in_tree.back().assign_branch_number_to_internal_nodes();
	// Print tree
	if (print_tree) {
	    if (read_trees >= tree_interval_start && read_trees <= tree_interval_end) {
		for (vector<tree>::iterator i=in_tree.begin(); i!=in_tree.end(); ++i) {
		    if (output_format == 'w') i->print_newick( print_br_length );
		    else if (output_format == 'x') {
			if (i==in_tree.begin() && ((tree_interval_start == 0 && read_trees == 1) || read_trees == tree_interval_start)) i->print_nexus_tree_intro(taxa_trans);
			stringstream ss;
			ss << "tree" << read_trees;
			if (i != in_tree.begin()) ss << '_' << (i - in_tree.begin());
			string tree_name(ss.str());
			i->print_tree_to_nexus( tree_name, print_br_length, true, taxa_trans );
		    }
		    else if (output_format == 's' && (svg_width < -0.5 || svg_height < -0.5)) i->print_svg_autoscale(true, true, fontsize, tip_colors);
		    else if (output_format == 's') i->print_svg(print_scalebar,print_nodelabel,svg_width,svg_height,textoffset,strokewidth,fontsize,font,tip_colors);
		}
	    }
	    if (read_trees > tree_interval_end) break;
	}
    }
    // Close file stream if open
    if (input_file.is_open()) input_file.close();
    if (print_tree && output_format == 'x' && read_trees >= tree_interval_start) cout << "End;" << endl;
    if (method == 'R') {
	std::cout << matrix.size() << ' ' << matrix.begin()->second.size() << std::endl;
	map<string, vector<char> >::const_iterator i;
	unsigned int max_length=0;
	for (i = matrix.begin(); i!= matrix.end(); ++i) {
	    unsigned int length = i->first.length();
	    if (length > max_length) max_length = length;
	}
	for (i = matrix.begin(); i!= matrix.end(); ++i) {
	    std::cout << i->first;
	    for (unsigned int n=0; n < max_length+1 - i->first.length(); ++n) std::cout << ' ';
	    for (vector<char>::const_iterator j = i->second.begin(); j!= i->second.end(); ++j) {
		std::cout << *j;
	    }
	    std::cout << std::endl;
	}
    }
}

void help () {
    std::cout << "Treebender is a command line program for manipulating trees." << endl;
    std::cout << "The program take a tree in newick format as indata through standard in." << endl;
    std::cout << "(c) Martin Ryberg 2014." << endl << endl;
    std::cout << "Usage:" << endl << "treebender [arguments] < file.tree" << endl;
    std::cout << "treebender [arguments] file.tree" << endl << endl;
    std::cout << "For the second alternative you need to be careful so treebender does not interpret the filename as an extra argument to a switch. If this happen treebender will expect" <<endl;
    std::cout << "input from standard in and it will appear as nothing is happening. This can be avoided by giving the filename after the switch --file/-f (see below)." <<endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--aMPL [text_string]                                       time calibrate tree with the adjusted Mean Path Length method (Svennblad 2008, Syst. Bio. 57:947-954). Give" << endl;
    std::cout << "                                                               calibration points with a text string with two comma separated taxa with the node as MRCA and the age" << endl;
    std::cout << "                                                               after colon. Several calibration points can be given separated by semicolon. The root node can be given" << endl;
    std::cout << "                                                               as root. eg. --aMPL \"taxon1,taxon2:50;root:100;taxon4,taxon5:20\" (default root:1)" << endl;
    std::cout << "--branch_lengths / -a [r/n]                                print branch lengths, the separator can be given as first argument after the switch, e.g. -a '\\n' (default" <<endl;
    std::cout << "                                                               is ','). If the switch r is given as second argument after the switch, e.g. -a '\\n' r, the value of the" << endl;
    std::cout << "                                                               root branch will be printed as well, if n (default) is given it will not." << endl;
    std::cout << "--change_names / -c [taxon_string]                         change the name of tips. Tip names should be pairwise with present name first and new name after separated" << endl;
    std::cout << "                                                               by '|'. Separate pairs should be separated by ',', e.g. -c 'taxon1|new1,taxon2|new2' (quotation marks" << endl;
    std::cout << "                                                               required). If several tips have the same name all will be changed. Changes later in the list will be" << endl;
    std::cout << "                                                               effected by changes made earlier in the list, e.g. -c 'taxon1|new1,new1|new2' will change the name of" << endl;
    std::cout << "                                                               taxon1 to new2." << endl;
    std::cout << "--clocklikness                                             gives the Z value that the two clades descending from each node have the same molecular clock." << endl;
    std::cout << "--depth / -D                                               get the distance from the root to the tip furthest away." << endl;
    std::cout << "--distances_to_root / -z [value_separator row_separator]   output the number of nodes and branch length distance to the root for each taxa. The separator between taxa" << endl;
    std::cout << "                                                               name and value can be specified, e.g. -p , \" | \" (dafault is newline and tab)." << endl;
    std::cout << "--drop_tips / -d [taxon_string]                            drop the given tips from the tree, e.g. -d taxon1,taxon2,taxon3" << endl;
    std::cout << "--get_tip_names / -t [separator]                           get the names of the tips/taxa in the tree, a separator can be specified, e.g. -t \\n (each name on" << endl;
    std::cout << "                                                               separate rows). ',' is the default separator." << endl;
    std::cout << "--get_branch_numbers                                       assign branch number as node label." << endl;
    std::cout << "--get_relaxed_outgroup [taxon_string]                      get the taxa in the clade that include the largest fraction of the difference between number of taxa included in" << endl;
    std::cout << "                                                               the given group and number not included in the group divided by the total number in the group." << endl;
    std::cout << "--file/-f                                                  give file name, e.g. -f file.tree." << endl;
    std::cout << "--format [newick/nexus]                                    give format of input, e.g. --format nexus. If no format is given and the input is a file treebender will try to" << endl;
    std::cout << "                                                               guess the format, if it is through standard in it will assume newick format." << endl;
    std::cout << "--help / -h                                                print this help." << endl;
    std::cout << "--internal_node_stats [value/param]                        print stats about the values on the intenal nodes. Counts nodes with value above given value, e.g." << endl;
    std::cout << "                                                               --internal_node_stats 1.96 (default: 1.0). If extra stats are given in FigTree/treeanotator format" << endl;
    std::cout << "                                                               the paramete to summarize can be given behind slash, i.e. --internal_node_stats 1.96/rate or" << endl;
    std::cout << "                                                               --internal_node_stats /rate" << std::endl;
    std::cout << "--interval [integer,integer]                               only print trees in the interval. Intervall given as first tree to print - last tree to print, e.g. --interval" << endl;
    std::cout << "                                                               10-100, or just the first tree to print." << endl;
    std::cout << "--inverse / -i                                             inverse the taxon string, e.g. drop all tips but the given" << endl;
    std::cout << "--is_monophyletic                                          test if the given taxa form a monophyletic group, e.g. --is_monophyletic taxon1,taxon2,taxon3" << endl;
    std::cout << "--ladderize / -l                                           laddrize the tree. If followed by l left ladderize, if followed by r right ladderize (default), e.g. -l r." << endl;
    std::cout << "--matrix_representation                                    present a phylip-formated matrix with splits of trees coded as characters. Intended for matrix representation" << endl;
    std::cout << "                                                               parsimony." << endl;
    std::cout << "--mid_point_root / -m                                      root the tree at the mid point." << endl;
    std::cout << "--multiply_branch_lengths / -u [value]                     multiply each branch in the tree with the given value, e.g. 3.5 (default 1.0)." << endl;
    std::cout << "--multiply_branch_lengths_clade / -V [value,taxon_string]  multiply branches in clades defined by the most recent common ancestor of comma separated taxa. Separate clades" << endl;
    std::cout << "                                                               with ':' or ';'. E.g. 3:Taxon_1,Taxon_2:Taxon_3,Taxon_4." << endl;
    std::cout << "--multiply_branch_lengths_until / -U [cut off,value]       multiply branches in tree up until cut off value distance from root with given value, e.g. 40,2 (default 0.0,1.0)." << endl;
    std::cout << "--n_supported [value]                                      get the number of nodes with higher support than given value. Should be followed by value, e.g. --n_supported 70.0" << endl;
    std::cout << "--newick / -w                                              output tree in newick format (default)." << endl;
    std::cout << "--nexus / -x                                               output tree in nexus format." << endl;
    std::cout << "--nni [integer/all]                                        perform nearest neighbor interchange. If a integer is given as extra argument the interchange will be done on that" << endl;
    std::cout << "                                                               branch (use --get_branch_numbers to get branch numbers). If 0 or no extra argument is given a branch will be" << std::endl;
    std::cout << "                                                               selected. If 'all' is given NNI will be performed for all branches, e.g. --nni 4 or --nni all." << std::endl;
    std::cout << "--no_branch_length / -0                                    do not print branch lengths. If there are no branch lengths in input tree the default is to print zero length" << endl;
    std::cout << "                                                               branches in the out tree. This argument override this and print no branch lengths." << endl;
    std::cout << "--number_of_taxa / -n                                      get the number of tips/taxa in the tree." << endl;
    std::cout << "--outgroup_root / -o [taxon_string]                        root using most recent common ancestor of given taxa, e.g. -o taxa1,taxa2." << endl;
    std::cout << "--patristic_distances / -p [value_separator row_separator] get the total patristic distance to all other taxa in the tree for each taxon, the separator between different" << endl;
    std::cout << "                                                               taxa, and the separator between taxon name and value can be specified, e.g. -p , \" | \" (default is new" << endl;
    std::cout << "                                                               line and tab)." << endl;
    std::cout << "--random_tree / -r                                         get a random topology (no branch lengths) with given number of taxa, e.g. -r 20 (default 0)." << endl;
    std::cout << "--read_figtree_annotations                                 will read annotations in FigTree/treeanotator format (e.g. [&rate=1.0,height=3.0])" << endl;
    std::cout << "--relaxed_outgroup_root [taxon_string]                     will root on the group defined as for --get_relaxed_outgroup." << endl;
    std::cout << "--set_branch_lengths / -b [value]                          set all branches in the tree to the given value, e.g. 0.5 (default 1.0)." << endl;
    std::cout << "--sum_branch_length / -s                                   get the sum of the branch lengths in the tree (including root branch if length for this is given)." << endl;
    std::cout << "--svg / -g                                                 output tree as svg immage. Extra graphical commands can be given as next argument. Ech command should be separated" << endl;
    std::cout << "                                                               by & and commands and arguments should be separated by :. Possible commands are: 'width' set width of figure;" << endl;
    std::cout << "                                                               'height' set hight of figure; 'offset' set offset between tips and tip label; 'strokewidth' set the width of" << endl;
    std::cout << "                                                               the branches; 'fontsize' sets the size of the font used; 'font' sets which font to use; and 'tipcolor' sets the" << endl;
    std::cout << "                                                               of the tip labels given in parenthesis directly behind the color. 'width' and 'height' are mandatory if you want to set" << endl;
    std::cout << "                                                               any other parameter than tip color. E.g. --svg 'width:300&height:400&tipcolor:red(taxon1,taxon2,taxon3)green(taxon4)'." << endl;
    std::cout << endl;

}

string pars_ARGV_string ( char input [] ) {
    string argv_string;
    bool file_pars(false);
    for (unsigned int i(0); input[i] != '\0'; ++i) {
	if (input[i] == ':' && !argv_string.compare("file")) {
	    file_pars = true;
	    argv_string.clear();
	}
	else argv_string += input[i];
    }
    if (file_pars) {
	ifstream fileinput;
	fileinput.open(argv_string.c_str());
	argv_string.clear();
	if (fileinput.is_open()) {
	    char in;
	    while (fileinput.get(in)) {
		//fileinput.get(in);// >> in;
		if (in != '\n' && in != '\r') argv_string += in;
	    }
	    if (fileinput.is_open()) fileinput.close();
	}
    }
    return argv_string;
}

