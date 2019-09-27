/********************************************************************
Copyright (C) 2016 Martin Ryberg

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
#include <vector>
#include <limits.h>
#include "tree.h"
#include "file_parser.h"
#include "clustertree.h"
#include "argv_parser.h"

using namespace std;

void help ();

int main (int argc, char *argv []) {
    #if DEBUG
    cerr << "Debugging version!!!" << endl;
    #endif //DEBUG
    char method('0');
    // from treesplitter
    const int max_trees = 1000;
    char split_criteria('n');
    char split_stop('s');
    unsigned int max_size(0);
    bool rooted(false);
    ///
    char tree_source('s');
    char output_format('w');
    char inverse_taxa('n');
    char randBr('y');
    unsigned int tree_interval_start(0);
    unsigned int tree_interval_end(UINT_MAX);
    string taxastring;
    vector<string> taxon_vector;
    vector< pair<float,float> > skyline_cut_offs;
    string separator(",");
    string separator2(" ");
    string tree_separator;
    bool print_scalebar(true);
    bool print_nodelabel(true);
    // svg
    float svg_width(-1.0);
    float svg_height(-1.0);
    float textoffset(5.0);
    int strokewidth(1);
    unsigned int fontsize(5);
    string font("Arial");
    string tip_colors;
    /////
    string stats_param;
    bool read_figtree_annotations = false;
    bool print_br_length(true);
    bool flag(false);
    bool right_ladderize(true);
    bool read_phylommand_block = false;
    float value(1.0);
    float cut_off(0.0);
    int n_taxa(0);
    unsigned int n_char(0);
    unsigned int n_rand_trees(1);
    unsigned int branch_no(0);
    map<string, vector<char> > matrix;
    string file_name;
    string file_format;
    ifstream input_file;
    file_parser input(&cin);
    /// from clustertree
    char cluster_method = '0';
    unsigned int name_position = 1;
    #ifdef DATABASE
    string database_data;
    #endif /*DATABASE*/
    bool quiet(true);
    /// Pars arguments ///
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--mid_point_root")) method = 'm';
            else if (!strcmp(argv[i],"-o") || !strcmp(argv[i],"--outgroup_root")) {
                method = 'g';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = argv_parser::pars_string( argv[++i]);
		else { cerr << "-o/--outgroup_root need an outgroup as next argument" << endl; return 1; }
            }
	    else if (!strcmp(argv[i],"--relaxed_outgroup_root")) {
                method = 'G';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = argv_parser::pars_string( argv[++i] );
                else { cerr << "--relaxed_outgroup_root need an outgroup as next argument" << endl; return 1; }
            }
            else if (!strcmp(argv[i],"--get_relaxed_outgroup")) {
                method = 'H';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = argv_parser::pars_string( argv[++i]);
		else { cerr << "--get_relaxed_outgroup need an outgroup as next argument" << endl; return 1; }
            }
            else if (!strcmp(argv[i],"--is_monophyletic")) {
                method = 'M';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = argv_parser::pars_string( argv[++i]); //argv[++i];
		else { cerr << "--is_monophyletic need a comma separated string of tip names to check for monophyly." << endl; return 1; }
            }
            else if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--ladderize")) {
                method = 'l';
                if ( i < argc-1 && argv[i+1][0] != '-') {
                    if (!(strcmp(argv[i+1],"l")) or !(strcmp(argv[i+1],"L")) or !(strcmp(argv[i+1],"left"))) { right_ladderize = 0; ++i; }
                    else if (!(strcmp(argv[i+1],"r")) or !(strcmp(argv[i+1],"R")) or !(strcmp(argv[i+1],"right"))) ++i;
                }
            }
	    else if (!strcmp(argv[i],"--split")) { // implement functionality previously in treesplitter
                if ( i < argc-1 && argv[i+1][0] != '-') {
////////////////////////////////////////////
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments);
		    if (!arguments[0].compare("l") || !arguments[0].compare("longest_branch")) {
			method = 'L';
			if (arguments.size() > 1 && !arguments[1].empty()) {
			    if ( arguments[1][0] == 'n' ) split_criteria = arguments[1][0];
			    else if ( arguments[1][0] == 'l' ) split_criteria = 'l';
			    else {
				cerr << "Unrecognized criterion for picking tree to split " << arguments[1] << "." << endl;
				return 1;
			    }
			}
			else split_criteria = 'l';
		    }
		    else if (!arguments[0].compare("m") || !arguments[0].compare("mid_point")) {
			method = 'P';
			if (arguments.size() > 1 && !arguments[1].empty()) {
			    if ( arguments[1][0] == 'n' ) split_criteria = arguments[1][0];
			    else if ( arguments[1][0] == 'l' ) split_criteria = 'l';
			    else {
				cerr << "Unrecognized criterion for picking tree to split " << arguments[1] << "." << endl;
				return 1;
			    }
			}
			else split_criteria = 'n';
		    }
		    else {
			cerr << "Method " << arguments[0] << " unrecognized." << endl;
			return 1;
		    }
///////////////////////////////////////
		}
		else {
		    cerr << "--split require a method for splitting the tree, see help for further instructions (--help)." << endl;
		}
	    }
	    else if (!strcmp(argv[i],"--split_stop")) { // implement functionality previously in treesplitter
                if ( i < argc-1 && argv[i+1][0] != '-') {
////////////////////////////////////////////
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments);
		    if (!arguments[0].compare("t") || !arguments[0].compare("max_tree_number")) {
			split_stop = 't';
		    }
		    else if (!arguments[0].compare("s") || !arguments[0].compare("max_tree_size")) {
			split_stop = 's';
		    }
		    else {
			cerr << "Unrecognized stop criterion: " << arguments[0] << "." << endl;
		    }
		    if ( arguments.size() > 1 && !arguments[1].empty() ) max_size = atoi(arguments[1].c_str());
		}
		else {
		    cerr << "--split_stop require a criterion for stopping splitting the tree, see help for further instructions (--help)." << endl;
		}
	    }
	    ///////////////a
	    //////// Implement functionality previously in clustertree
	    else if (!strcmp(argv[i],"--cluster")) {
		method = 'C';
		if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
		    if (!arguments[0].compare("long_branch")) {
			cluster_method = 'l';
			if (arguments.size() > 1) cut_off = atof(arguments[1].c_str());
		    }
		    else if (!arguments[0].compare("tip_name")) {
			cluster_method = 't';
	    		if (arguments.size() > 1 && !arguments[1].empty()) {
			    separator = arguments[1][0];
			    if (arguments.size() > 2 && !arguments[2].empty()) {
				name_position = atoi(arguments[2].c_str());
			    }
			}
		    }
		    else if (!arguments[0].compare("branch_length")) {
			cluster_method = 'b';
			if (arguments.size() > 1) cut_off = atof(arguments[1].c_str());
		    }
		    #ifdef DATABASE
		    else if (!arguments[0].compare("database")) {
			cluster_method = 'd';
			if (arguments.size()>1 && !arguments[1].empty()) database_data = argv[++i];
			else {
			    cerr << "database clustering require a database file name as next argument." << endl;
			    return 1;
			}
		    }
		    #endif //DATABASE
		    else {
			cerr << "Unknown clustering method: " << arguments[0] << endl;
			return 1;
		    }
		}
		else {
		    cerr << "--cluster require a clustering method as next argument." << endl;
		}
	    }
	    else if (!strcmp(argv[i],"--cut_off")) cut_off = atof(argv[++i]);
    	    else if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose")) quiet = false;
	    else if (!strcmp(argv[i],"-R") || !strcmp(argv[i],"--rooted")) rooted=true;
            else if (!strcmp(argv[i],"-a") || !strcmp(argv[i],"--branch_lengths")) {
                method = 'a';
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
                    separator = arguments[0];
		    if (arguments.size()>1) {
                        if (arguments[1][0] == 'r') { flag=true; ++i; }
                        else if (arguments[1][0] == 'n') { flag=false; ++i; }
		    }
		    if (arguments.size()>2) {
			tree_separator = arguments[2];
		    }
/*
                    separator = argv[++i];
                    if (i < argc-1 && argv[i+1][0] != '-') {
                        if (argv[i+1][0] == 'r' && strlen(argv[i+1])==1) { flag=true; ++i; }
                        else if (argv[i+1][0] == 'n' && strlen(argv[i+1])==1) { flag=false; ++i; }
                    }*/
                }
            }
	    else if (!strcmp(argv[i],"--clade_credibility"))  method = '*';
	    else if (!strcmp(argv[i],"--clear_internal_node_labels")) method = '-';
            else if (!strcmp(argv[i],"-z") || !strcmp(argv[i],"--distances_to_root")) {
                method = 'o';
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
		    separator = arguments[0];
		    if (arguments.size() > 1 && !arguments[1].empty()) {
	    		separator2 = arguments[1];
		    }
                    else separator2 = "\t";
		    if (arguments.size() > 2) {
	    		tree_separator = arguments[2];
		    }
                }
                else separator = "\n";
            }
            else if (!strcmp(argv[i],"-p") || !strcmp(argv[i],"--patristic_distances")) {
                method = 'p';
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
		    separator = arguments[0];
		    if (arguments.size() > 1 && !arguments[1].empty()) {
	    		separator2 = arguments[1];
		    }
                    else separator2 = "\t";
		    if (arguments.size() > 2) {
	    		tree_separator = arguments[2];
		    }
                    /*separator = argv[++i];
                    if ( i < argc-1 && argv[i+1][0] != '-') {
                        separator2 = argv[++i];
                    }*/
                }
                else separator = "\n";
            }
            else if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"--set_branch_lengths")) {
                method = 'b';
                if ( i < argc-1 && argv[i+1][0] != '-') value = atof(argv[++i]);
            }
	    else if (!strcmp(argv[i],"--null_short_branches")) {
		method = '<';
		if ( i < argc-1 && argv[i+1][0] != '-') {
		    value = atof(argv[++i]);
		}
		else {
		    cerr << "--null_short_branches require a value for what should be considered a short branch as next argument";
		}
	    }
            else if (!strcmp(argv[i],"-u") || !strcmp(argv[i],"--multiply_branch_lengths")) {
                method = 'u';
                if ( i < argc-1 && argv[i+1][0] != '-') value = atof(argv[++i]);
            }
	    else if (!strcmp(argv[i],"-U") || !strcmp(argv[i],"--multiply_branch_lengths_until")) {
		method = 'U';
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
		    value = atof(arguments[0].c_str());
		    if (arguments.size() > 1 && !arguments[1].empty()) {
	    		cut_off = atof(arguments[1].c_str());
		    }
		}
	    }
	    else if (!strcmp(argv[i],"-S") || !strcmp(argv[i],"--multiply_branch_lengths_skyline")) {
		method='$';
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ',', arguments );
		    for (vector<string>::iterator i=arguments.begin(); i != arguments.end(); ++i) {
			vector<string> second_level;
			argv_parser::pars_sub_args(i->c_str(), ':', second_level);
			skyline_cut_offs.push_back(make_pair(atof(second_level[0].c_str()),atof(second_level[1].c_str())));
		    }
		    #ifdef DEBUG
	    	    cerr << "Number of cut offs " << skyline_cut_offs.size() << endl;
    		    #endif //DEBUG
		}
		else { cerr << "-S/--multiply_branch_lengths_skyline require a list of time cut offs and rate multipliers as second argument." << endl; return 1; }
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
                if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
		    separator = arguments[0];
		    if (arguments.size() > 1) {
	    		tree_separator = arguments[1];
		    }
		}
                //if ( i < argc-1 && argv[i+1][0] != '-') separator = argv[++i];
            }
            else if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--change_names")) {
                method = 'c';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = argv_parser::pars_string( argv[++i]);
            }
            else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--drop_tips")) {
                method = 'd';
                if ( i < argc-1 && argv[i+1][0] != '-' ) taxastring = argv_parser::pars_string( argv[++i]);
            }
	    else if (!strcmp(argv[i],"--drop_tips_on_short_branches")) {
		method = '>';
                if ( i < argc-1 && argv[i+1][0] != '-' ) cut_off = atof(argv[++i]);
		else { cerr << "--drop_tips_on_short_branches require a number to use as cut off as next argument." << endl; return 1; } 
	    }
            else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--random_tree")) {
                tree_source = 'r';
                if ( i < argc-1 && (argv[i+1][0] == '0' || argv[i+1][0] == '1' || argv[i+1][0] == '2'
		   || argv[i+1][0] == '3' || argv[i+1][0] == '4' || argv[i+1][0] == '5'
		   || argv[i+1][0] == '6' || argv[i+1][0] == '7' || argv[i+1][0] == '8'
		   || argv[i+1][0] == '9' || argv[i+1][0] == ':') ) {
			++i;
			vector<string> arguments;
			argv_parser::pars_sub_args(argv[i], ':', arguments );
			if (!arguments[0].empty()) n_taxa = atoi(arguments[0].c_str());
			if (arguments.size() > 1 && !arguments[1].empty()) {
			    n_rand_trees = atoi(arguments[1].c_str());
			}
			//n_taxa = atoi(argv[++i]);
		    }
            }
	    else if (!strcmp(argv[i],"--rand_node_depths")) { // || !strcmp(argv[i],"--yule")) {
		method = 'B';
		if (i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    if (argv[i][0] == 'Y' || argv[i][0] == 'y')	randBr='y';
		    else if (argv[i][0] == 'C' || argv[i][0] == 'c') randBr='c';
		}
	    }
	    else if (!strcmp(argv[i],"--gamma")) method = '!';
            else if (!strcmp(argv[i],"--interval")) {
                if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    string argument(argv[++i]);
		    size_t separator = argument.find("-");
		    tree_interval_start = atoi(argument.substr(0,separator).c_str());
		    if (separator < argument.length()) tree_interval_end = atoi(argument.substr(separator+1).c_str());
		    if (tree_interval_end < tree_interval_start) {
			cerr << "The end of the interval (--interval) has to be after the start of the interval." << endl;
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
                method = 'A';
                if ( i < argc-1 && argv[i+1][0] != '-') taxastring = argv_parser::pars_string( argv[++i] );
            }
            else if (!strcmp(argv[i],"--nni")) {
                method = 'I';
		if (i+1 < argc && argv[i+1][0] != '-') {
		    ++i;
		    if (!strcmp(argv[i],"all")) flag = true;
		    else branch_no = atoi(argv[i]);
		}
	    }
            else if (!strcmp(argv[i],"--clocklikness")) {
		method = 'K';
	    }
	    else if (!strcmp(argv[i],"--n_char")) {
                if ( i < argc-1 && argv[i+1][0] != '-') n_char = atoi(argv[++i]);
		else { cerr << "--n_char require the number of characters as next argument." << endl; return 1; }
	    }
	    else if (!strcmp(argv[i],"--read_figtree_annotations")) read_figtree_annotations = true;
            else if (!strcmp(argv[i],"--matrix_representation")) method = 'R';
            else if (!strcmp(argv[i],"--internal_node_stats")) {
		method = 'T';
		if ( i < argc-1 && argv[i+1][0] != '-') {
		    ++i;
		    vector<string> arguments;
		    argv_parser::pars_sub_args(argv[i], ':', arguments );
		    if (!arguments[0].empty()) value = atof(arguments[0].c_str());
		    if (arguments.size() > 1 && !arguments[1].empty()) {
	    		stats_param = arguments[1];
			read_figtree_annotations = true;
		    }
                }
		/*    ++i;
		    int j=0;
		    while (argv[i][j] != '\0' && argv[i][j] != '/') ++j;
		    string temp;
		    for (int k=0; k<j; ++k) temp+=argv[i][k];
		    value = atof(temp.c_str());
		    if (argv[i][j] != '\0') ++j;
		    while (argv[i][j] != '\0') stats_param += argv[i][j++];
		}*/
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
            else if (!strcmp(argv[i],"--output")) {
		if ( i < argc-1 && argv[i+1][0] != '-' ) {
		    ++i;
		    if (!strcmp(argv[i],"nexus") || !strcmp(argv[i],"nex") || (argv[i][0] == 'x' && argv[i][1] == '\0')) output_format = 'x';
		    else if (!strcmp(argv[i],"newick") || !strcmp(argv[i],"new") || (argv[i][0] == 'w' && argv[i][1] == '\0')) output_format = 'w';
		    else if ((argv[i][0] == 's' || argv[i][0] == 'S') && (argv[i][1] == 'v' || argv[i][1] == 'V') && (argv[i][2] == 'g' || argv[i][2] == 'G') && (argv[i][3] == '\0' || argv[i][3] == ':')) {
			output_format = 's';
			unsigned int j=3;
		    	string temp;
	    		string param;
    			while (1) {
			    if ((argv[i][j] == '&' || argv[i][j] == '\0') && !param.empty()) {
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
			    else if (argv[i][j] == '\0') break;
			    else if (argv[i][j] == ':') {
				param = temp;
				temp.clear();
			    }
			    else temp += argv[i][j];
			    ++j;
			}
		    }
		    else {
			std::cerr << "Do not recognize format " << argv[i] << "." << endl;
			return 1;
		    }
		}
		else std::cerr << "--output require nexus(nex or x), newick (new or w), or svg (SVG) as additional argument" << endl;
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
    if (!quiet) {
       std::cerr << "The program was called with the following command:" << endl;
       for (int i=0; i<argc; ++i) std::cerr << argv[i] << ' ';
       std::cerr << endl << endl;
    }
    /// Input error check
    if ( (method == 'c' || method == 'd' || method == 'M') && taxastring.empty() ) {
        std::cerr << "Selected method reqire a taxon string. Use -h to print help." << endl;
        return 1;
    }
    #ifdef DEBUG
    if (!taxastring.empty()) {
	cout << "Taxon string: " << taxastring << endl;
    }
    #endif //DEBUG
    /// Enterprete escape
    if ( !separator.compare("\\n") ) separator = "\n";
    else if ( !separator.compare("\\t") ) separator = "\t";
    else if ( !separator.compare("\\r") ) separator = "\r";
    /// prepairing to read trees
    if (!file_name.empty()) {
	#ifdef DEBUG
	cerr << "File name " << file_name << endl;
	#endif //DEBUG
	input_file.open(file_name.c_str(),std::ifstream::in);
	if (input_file.good())
	    input.file_stream = &input_file;
	else {
	    cerr << "Could not open file: " << file_name << endl;
	    return 1;
	}
	if (!quiet) cerr << "Taking input from " << file_name << "." << endl;
    }
    /// Setting variables
    bool print_tree = true;
    map <string,string> taxa_trans;
    #ifdef DEBUG
    cerr << "File type " << input.get_file_type() << endl;
    #endif //DEBUG
    /// Setting file format
    if (!file_format.empty()) {
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
    if (!quiet) cerr << "File format: " << input.get_file_type() << endl;
    /// If doing aMPL and nexus is given check data in phylommand block
    if (method=='A' && !taxastring.compare("nexus")) {
	read_phylommand_block = true;
	taxastring.clear();
    }
    if (read_phylommand_block) {
	if (!quiet) cerr << "Reading phylommand block" << endl;
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
		if (method=='A') {
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
    // Prepare to read trees if nexus
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
    unsigned int read_trees(0);
////// Cluster tree ///////
    #ifdef DATABASE
    string database;
    string table;
    string key;
    string column;
    if (cluster_method == 'd' && database_data.empty()) {
	if (!quiet) {
	    cerr << "If clustering using database entries the database, table, column for tree tip" << endl;
	    cerr << "annotation, and column used for clustering must be given as a comma separated" << endl;
	    cerr << "string, e.g. database,table,accno_column,species_column" << endl;
	}
	return 1;
    }
    else if (cluster_method == 'd') {
	unsigned int k=0;
	for (unsigned int i=0; i < database_data.length(); ++i) {
	    if (database_data[i] == ',') ++k;
	    else if (k == 0) database += database_data[i];
	    else if (k == 1) table += database_data[i];
	    else if (k == 2) key += database_data[i];
	    else if (k == 3) column += database_data[i];
	}
    }
    #endif /*DATABASE*/
    if (cut_off < 0) {
	if (!quiet) cerr << argv[0] << "Only positive cut off values are accepted." << endl;
	return 1;
    }
///////////////////////////////////////
    /// Starting to read trees
    if (!quiet) cerr << "Starting to process trees" << endl;
    while (1) {
	vector<tree> in_tree;
	in_tree.push_back(tree());
	// Get tree
	if (tree_source == 's') {
	    if (input.test_file_type("nexus") || input.test_file_type("newick")) {
		if (input.test_file_type("nexus")) {
		    if (read_trees != 0) nexus_command = input.read_next_nexus_command();
		    if (!(nexus_command==nexus_command::TREE && input.move_to_start_of_tree()))
			break;
		}
		if (method == 'C') { // clustertree
		    clustertree tree;
		    tree.tree_file_parser( *(input.file_stream), taxa_trans, read_figtree_annotations );
		    ++read_trees;
    		    if (tree.empty()) break;
		    cout << "### tree " << read_trees << " ###" << endl;
		    if (cluster_method == 'l') tree.br_length_clust_max_cout( cut_off );
		    else if (cluster_method == 'b') tree.short_br_clust ( cut_off );
		    else if (cluster_method == 't') tree.name_clust_cout( name_position, separator[0]);
		    #ifdef DATABASE
		    else if (cluster_method == 'd') {
			cout << "[" << database << ", " << table << ", " << key << ", " << column << "]" << endl;
			tree.db_clust_cout( database.c_str(), table, key, column);
		    }
		    #endif /*DATABASE*/
		    continue;
		}
		else {
		    in_tree.back().tree_file_parser( *(input.file_stream), taxa_trans, read_figtree_annotations );
		}
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
	if (in_tree.back().empty()) break;
	if (read_trees < tree_interval_start) continue;
	if (read_trees > tree_interval_end) break;
	// Make preparations
	if (inverse_taxa == 'y') taxastring = in_tree.back().not_present( taxastring );
    
	// Do manipulation
	if (method == 'm') in_tree.back().midpoint_root( );
	else if (method == 'M') {
	    if (!quiet) cerr << "Is monophyletic?:" << endl;
	    if ( in_tree.back().is_monophyletic( taxastring ) ) cout << "Yes" << endl;
	    else cout << "No" << endl;
	    print_tree = false;
	}
	else if (method == 't') {
	    cout << tree_separator;
	    if (!quiet) cerr << "Tip names:" << endl;
	    cout << in_tree.back().tip_names( separator ) << endl;
	    print_tree = false;
	}
	else if (method == '-') in_tree.back().clear_internal_node_labels();
	else if (method == 'g') in_tree.back().outgroup_root( taxastring );
	else if (method == 'G') in_tree.back().relaxed_outgroup_root( taxastring );
	else if (method == 'H') {
	    if (!quiet) cerr << "Outgroup:" << endl;
	    in_tree.back().print_relaxed_outgroup( taxastring );
	    cout << endl;
	    print_tree = false;
	}
	else if (method == 'o') {
	    cout << tree_separator;
	    if (!quiet) cerr << "Distances:" << endl;
	    in_tree.back().print_distance_to_root_for_all ( separator2, separator );
	    print_tree = false;
	}
	else if (method == 'l') in_tree.back().ladderize(right_ladderize);
	else if (method == 'd') in_tree.back().drop_tips( taxastring );
	else if (method == '>') in_tree.back().drop_tips_on_short_branches (cut_off);
	else if (method == 'c') in_tree.back().change_tip_names( taxastring );
	else if (method == 'n') {
	    if (!quiet) cerr << "Number of tips:" << endl;
	    cout << in_tree.back().n_tips() << endl;
	    print_tree = false;
	}
	else if (method == 'D') {
	    //cout << "### tree " << read_trees << " ###" << endl;
	    if (!quiet) cerr << "Longest distance between root and tip:" << endl;
	    cout << in_tree.back().longest_to_tip() << endl;
            print_tree = false;
	}
	else if (method == 's') {
	    if (!quiet) cerr << "Tree length:" << endl;
	    cout << in_tree.back().sum_br_length() << endl;
	    print_tree = false;
	}
	else if (method == '*') {
	    if (!quiet) cerr << "Log product of support values:" << endl;
	    cout << in_tree.back().log_clade_credibility() << endl;
	    print_tree = false;
	}
	else if (method == 'u') in_tree.back().multiply_br_length( value );
	else if (method == 'U') in_tree.back().multiply_br_length_cut_off( cut_off, value );
	else if (method == '$') in_tree.back().multiply_br_length_skyline(skyline_cut_offs);
	else if (method == 'V') in_tree.back().multiply_br_length_clades( taxon_vector, value );
	else if (method == '<') in_tree.back().null_short_branches(value);
	else if (method == 'b') in_tree.back().set_br_length( value );
	else if (method == 'a') {
	    cout << tree_separator;
	    if (!quiet) cerr << "Branch lengths:" << endl;
	    in_tree.back().print_branch_lengths ( separator, flag );
	    print_tree = false;
	}
	else if (method == 'I') {
	    unsigned int start = 2; // not 1, only do branch swapping along the root branch once
	    unsigned int n_tips = in_tree.back().n_tips();
	    if (n_tips < 4) std::cerr << "Warning!!! Can not do NNI on tree with fewer than 4 tips." << std::endl;
	    else {
		if (branch_no > (n_tips-3))  std::cerr << "Warning!!! Branch to swap on is out of bound (" << branch_no << " > " << (n_tips-3) << ")." << std::endl;
		else {
		    if (flag) branch_no = (n_tips-3); // if swapping on all
		    else start = branch_no; // if only swapping on one branch
		    for (unsigned int i=start; i<=branch_no; ++i) {
			tree swapOne;
			tree swapTwo;
			in_tree[0].nni(i,swapOne,swapTwo); // if i == 0 choose random tip
			in_tree.push_back(swapOne); // swap_one;
			in_tree.push_back(swapTwo); // swap_one;
		    }
		}
	    }
	}
	else if (method == 'p') {
	    cout << tree_separator;
	    if (!quiet) cerr << "Patristic distances: " << endl;
	    in_tree.back().print_distances_to_all_for_each_taxa ( separator2, separator );
	    print_tree = false;
	}
	else if (method == 'Y') {
	    if (!quiet) cerr << "Number of nodes with support above " << value << ":" << endl;
	    cout << in_tree.back().n_supported(value) << endl;
	    print_tree = false;
	}
	else if (method == '!') {
	    if (!quiet) cerr << "Gamma: " << endl;
	    cout << in_tree.back().gamma() << endl;
	    print_tree = false;
	}
       	else if (method == 'A') {
	    if (taxastring.empty()) taxastring="root:1";
	    if (n_char < 1) n_char = 1;
	    #ifdef DEBUG
	    cerr << taxastring << endl;
	    #endif //DEBUG
	    in_tree.back().adjustedMPL(taxastring, n_char);
	}
       	else if (method == 'K') {
	    in_tree.back().test_clock_likness(n_char);
	}
	else if (method == 'T') {
	    cout << "### tree " << read_trees << " ###" << endl;
	    if (!quiet) cerr << "Statistics for " << stats_param << "(cut off " << value << "):" << endl;
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
	else if (method == 'L' || method == 'P') {
	    if (!quiet) cerr << "Starting to split tree " << read_trees << endl;
	    unsigned int split_count = 0;
	    vector<tree>::iterator tree_to_split = in_tree.begin();
	    bool cont = true;
	    if (max_size==0) max_size = tree_to_split->n_tips() + 100;
	    while(cont && in_tree.size() <= max_trees) {
		if (!quiet) cerr << "Making split " << ++split_count << endl;
		cont=false;
		tree tree_part;
		if (method == 'L') {
		    if (rooted) tree_to_split->split_tree_at_longest_branch(&tree_part);
		    else tree_to_split->split_tree_at_longest_branch_unrooted(&tree_part);
		}
		else if (method == 'P') tree_to_split->split_tree_at_midpoint(&tree_part);
		in_tree.push_back(tree_part);
		#ifdef DEBUG
		cerr << "Made split " << split_count << endl;
		#endif //DEBUG
		unsigned int max=0;
		double max_br=0;
		for (vector<tree>::iterator i=in_tree.begin(); i != in_tree.end(); ++i) {
		    unsigned int n = i->n_tips();
		    if (split_stop == 's')
			if (n > max_size) cont = true;
		    if (split_criteria == 'n') {
			if (n > max) {
			    max=n;
			    tree_to_split = i;
			}
		    }
		    else if (split_criteria == 'l') {
			double br=i->longest_branch();
			if (br > max_br) {
			    max_br=br;
			    tree_to_split = i;
			}
		    }
		}
		if (split_stop == 't')
		    if ( in_tree.size() < max_size ) cont = true;
	    }
	}
	else if (method == 'B') {
	    in_tree.back().set_br_length( 0.0 );
	    if (randBr == 'c') in_tree.back().add_Coal_node_depths();
	    else in_tree.back().add_Yule_node_depths();
	}
	// Print tree
	if (print_tree) {
	    //if (read_trees >= tree_interval_start && read_trees <= tree_interval_end) {
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
		    else if (output_format == 's' && (svg_width < -0.5 || svg_height < -0.5)) {
			//unsigned int n_tips = i->n_tips();
			//if (100/n_tips > fontsize) fontsize = 100/n_tips;
			i->print_svg_autoscale(true, true, fontsize, tip_colors);
		    }
		    else if (output_format == 's') i->print_svg(print_scalebar,print_nodelabel,svg_width,svg_height,textoffset,strokewidth,fontsize,font,tip_colors);
		}
	    //}
	    //if (read_trees > tree_interval_end) break;
	}
    }
    // Close file stream if open
    if (input_file.is_open()) input_file.close();
    if (print_tree && output_format == 'x' && read_trees >= tree_interval_start) cout << "End;" << endl;
    if (method == 'R') {
	map<string, vector<char> >::const_iterator i;
	for (i = matrix.begin(); i!= matrix.end(); ++i) {
	    cout << '>' << i->first << endl;
	    for (vector<char>::const_iterator j = i->second.begin(); j!= i->second.end(); ++j) {
		cout << *j;
	    }
	    cout << std::endl;
	}
    }
}

void help () {
    cout << "Treebender " << VERSION << " is a command line program for manipulating trees. The program" << endl;
    cout << "take a tree in newick or nexus format as indata through standard in or from a" << endl;
    cout << "file and process it according to given options." << endl;
    cout << "(c) Martin Ryberg " << YEAR << "." << endl << endl;
    cout << "Usage:" << endl << "treebender [arguments] < file.tree" << endl;
    cout << "treebender [arguments] file.tree" << endl << endl;
    cout << "For the second alternative you need to be careful so treebender does not" << endl;
    cout << "interpret the filename as an extra argument to a switch. If this happen" << endl;
    cout << "treebender will expect input from standard in and it will appear as nothing is" << endl;
    cout << "happening. This can be avoided by giving the filename after the switch --file/" << endl;
    cout << "-f. When taxa should be given as extra arguments they can be given in a text" << endl;
    cout << "following the format for the argument. Newline and carriage returns will be" << endl;
    cout << "ignored. The file name should be given behind the word file and colon, e.g." << endl;
    cout << "-d file:file_name.txt." <<endl << endl;
    cout << "Arguments:" << endl;
    /*cout << "--aMPL [string]                   time calibrate tree with the adjusted Mean" << endl;
    cout << "                                  Path Length method (Svennblad 2008, Syst. Bio." << endl;
    cout << "                                  57:947-954). Give  calibration points with a" << endl;
    cout << "                                  text string with two comma separated taxa with" << endl;
    cout << "                                  the node as MRCA and the age  after colon." << endl;
    cout << "                                  Several calibration points can be given" << endl;
    cout << "                                  separated by semicolon. The root node can be" << endl;
    cout << "                                  given as root. eg. --aMPL \"taxon1,taxon2:50;" << endl;
    cout << "                                  root:100;taxon4,taxon5:20\" (default root:1)." << endl; */
    cout << "--branch_lengths / -a [r/n]       print branch lengths, the separator can be" << endl;
    cout << "                                  given as first argument after the switch, e.g." <<endl;
    cout << "                                   -a '\\n' (default is ','). If the switch r is" << endl;
    cout << "                                  given as second argument after a colon (:)," << endl;
    cout << "                                  e.g. -a '\\n:r', the value of the root branch" << endl;
    cout << "                                  will be printed as well, if n (default) is" << endl;
    cout << "                                  given it will not.  A separator between output" << endl;
    cout << "                                  from different trees can be given after" << endl;
    cout << "                                  another colon." << endl;
    cout << "--change_names / -c [taxa]        change the name of tips. Tip names to be" << endl;
    cout << "                                  changed should be given pairwise with present" << endl;
    cout << "                                  name first and new name second, separated by" << endl;
    cout << "                                  '|'. Separate pairs should be separated by ','" << endl;
    cout << "                                  e.g. -c 'taxon1|new1,taxon2|new2' (quotation" << endl;
    cout << "                                  marks required). If several tips have the same" << endl;
    cout << "                                  name all will be changed. Changes later in the" << endl;
    cout << "                                  list will be effected by changes made earlier" << endl;
    cout << "                                  in the list, e.g. -c 'taxon1|new1,new1|new2'" << endl;
    cout << "                                  will change the name of taxon1 to new2." << endl;
    /*cout << "--clocklikness                    gives the Z value that the two clades" << endl;
    cout << "                                  descending from each node have the same" << endl;
    cout << "                                  molecular clock." << endl; */
    cout << "--clade_credibility               give the log of the product of the support of" << endl;
    cout << "                                  all clades." << endl;
    cout << "--clear_internal_node_labels      delete the internal node labels" << endl;
    cout << "--cluster [method]                get clusters based on method, e.g. --cluster" << endl;
    cout << "                                  branch_length. Available methods:" << endl;
    cout << "                                  branch_length - separate clusters by single" << endl;
    cout << "                                     link clustering based on phylogenetic" << endl;
    cout << "                                     distance. Cut off should be given after" << endl;
    cout << "                                     colon, e.g. --cluster branch_length:0.03." << endl;
    #ifdef DATABASE
    cout << "                                  database - cluster based on annotations" << endl;
    cout << "                                     available in SQLite database. Need to be" << endl;
    cout << "                                     followed by a comma separated string with the" << endl;
    cout << "                                     database, table, column for tree tip" << endl;
    cout << "                                     annotation, and column used for" << endl;
    cout << "                                     clustering given, e.g.: --cluster database:" << endl;
    cout << "                                     database_file,table,accno_column,"<< endl;
    cout << "                                     species_column." << endl;
    #endif /*DATABASE*/
    cout << "                                  long_branch - returns taxa in clades on" << endl;
    cout << "                                     branches longer than cut off. Cut off" << endl;
    cout << "                                     should be given after colon (:)." << endl;
    cout << "                                  tip_name - cluster taxa based on taxon" << endl;
    cout << "                                     annotation. Should be followed after a" << endl;
    cout << "                                     colon by a single character that" << endl;
    cout << "                                     separates different parts of the tip name" << endl;
    cout << "                                     (default ' ') and after another colon (:) a" << endl;
    cout << "                                     number giving which position in the name" << endl;
    cout << "                                     should be used for clustering, (default 1)," << endl;
    cout << "                                     e.g. tip_name:\\|:5." << endl;
    //cout << "--cut_off                         the cut off to use when clustering." << endl;
    cout << "--depth / -D                      get the longest distance from the root to any" << endl;
    cout << "                                  of the tips." << endl;
    cout << "--distances_to_root / -z [sep.]   output the number of nodes and branch length" << endl;
    cout << "                                  distance to the root for each tip. The" << endl;
    cout << "                                  separator between tip name and value can be" << endl;
    cout << "                                  specified, separated by colon, e.g. -p \",:|\"" << endl;
    cout << "                                  (default is newline and tab). A separator" << endl;
    cout << "                                  between output from different trees can be" << endl;
    cout << "                                  given after another colon" << endl;
    cout << "--drop_tips / -d [taxa]           drop the given tips from the tree, e.g. -d" << endl;
    cout << "                                  taxon1,taxon2,taxon3." << endl;
    cout << "--drop_tips_on_short_branches     drop tips on branches shorter than given cut" << endl;
    cout << "                                  off (given as next argument). Tips to the left" << endl;
    cout << "                                  are dropped before tips to the right." << endl;
    cout << "--gamma                           get the gamma value (Pybus & Harvey 2000," << endl;
    cout << "                                  Proc. R. Soc. Lond. B 267, 2267-2272)"<< endl;
    cout << "--get_tip_names / -t [sep.]       get the names of the tips in the tree, a" << endl;
    cout << "                                  separator can be specified, e.g. -t \\\\n (each" << endl;
    cout << "                                  name on separate rows; ',' is the default" << endl;
    cout << "                                  separator).  A separator between output from" << endl;
    cout << "                                  different trees can be given after another" << endl;
    cout << "                                  colon" << endl;
    cout << "--get_branch_numbers              assign branch numbers as node labels." << endl;
    cout << "--get_relaxed_outgroup [taxa]     get the taxa in the clade that include the" << endl;
    cout << "                                  largest fraction of the difference between" << endl;
    cout << "                                  number of taxa included in the given group and" << endl;
    cout << "                                  number not included in the group divided by" << endl;
    cout << "                                  the total number in the group. Taxa given as" << endl;
    cout << "                                  comma separated string (see --drop_tips)." << endl;
    cout << "--file / -f [file]                give file name, e.g. -f file.tree." << endl;
    cout << "--format [newick/nexus]           give format of input, e.g. --format nexus. If" << endl;
    cout << "                                  no format is given and the input is a file" << endl;
    cout << "                                  treebender will try to guess the format, if it" << endl;
    cout << "                                  is through standard in it will assume newick" << endl;
    cout << "                                  format." << endl;
    cout << "--help / -h                       print this help." << endl;
    cout << "--internal_node_stats [val./par.] print stats about the values on the internal" << endl;
    cout << "                                  nodes. Counts nodes with value above given" << endl;
    cout << "                                  value, e.g. --internal_node_stats 1.96" << endl;
    cout << "                                  (default: 1.0). If extra stats are given in" << endl;
    cout << "                                  FigTree/treeanotator format the parameter to" << endl;
    cout << "                                  summarize can be given behind colon, e.g." << endl;
    cout << "                                  --internal_node_stats 1.96:rate, or" << endl;
    cout << "                                  --internal_node_stats :rate." << std::endl;
    cout << "--interval [integer-integer]      only print the trees in the interval. Interval" << endl;
    cout << "                                  given as first tree to print - last tree to" << endl;
    cout << "                                  print, e.g. --interval 10-100, or just the" << endl;
    cout << "                                  first tree to print, e.g. --interval 1000." << endl;
    cout << "--inverse / -i                    inverse the string of taxa, e.g. drop all tips" << endl;
    cout << "                                  but the given. E.g -d taxon1,taxon2,taxon3 -i" << endl;
    cout << "--is_monophyletic [taxa]          test if the given taxa form a monophyletic" << endl;
    cout << "                                  group, e.g. --is_monophyletic taxon1,taxon2." << endl;
    cout << "--ladderize / -l                  laddrize the tree. If followed by l - left" << endl;
    cout << "                                  ladderize, if followed by r - right ladderize" << endl;
    cout << "                                  (default), e.g. -l r." << endl;
    cout << "--matrix_representation           present a fasta-formated matrix with splits" << endl;
    cout << "                                  of trees coded as characters. Intended for" << endl;
    cout << "                                  matrix representation parsimony." << endl;
    cout << "--mid_point_root / -m             root the tree at the mid point." << endl;
    cout << "--multiply_branch_lengths /       multiply each branch in the tree with the" << endl;
    cout << "   -u [value]                     given value, e.g. 3.5 (default 1.0)." << endl;
    cout << "--multiply_branch_lengths_clade / multiply branches in clades defined by the" << endl;
    cout << "   -V [value:taxon_string]        most recent common ancestor of comma separated" << endl;
    cout << "                                  taxa. Separate clade with colon E.g. -V 3:" << endl;
    cout << "                                  Taxon_1,Taxon_2:Taxon_3,Taxon_4." << endl;
    cout << "--multiply_branch_lengths_until / multiply branches in tree up until given" << endl;
    cout << "   -U [value:cut off]             distance (cut off) from root with the given" << endl;
    cout << "                                  value (separated by colon), e.g. 2:40 (default" << endl;
    cout << "                                  1.0:0.0)." << endl;
    cout << "--multiply_branch_lengths_skyline like --multiply_branch_lengths_until but it is" << endl;
    cout << "/ -S                              posible to give comma separated list of rates" << endl;
    cout << "                                  and distances, where each rate will be applied" << endl;
    cout << "                                  between the previous distance and given" << endl;
    cout << "                                  distance. e.g. -S 2:40,1:50,10:60." << endl;
/*  cout << "--n_char                          give number of characters, e.g. --n_char 1000." << endl;
    cout << "                                  For example for --clocklikness and --aMPL." << endl;
*/
    cout << "--n_supported [value]             get the number of nodes with higher support" << endl;
    cout << "                                  than given value. Should be followed by value," << endl;
    cout << "                                  e.g. --n_supported 70.0." << endl;
    cout << "--nni [node/all]                  perform nearest neighbor interchange. If a" << endl;
    cout << "                                  integer is given as extra argument the" << endl;
    cout << "                                  interchange will be done on that branch (use" << endl;
    cout << "                                  --get_branch_numbers to get branch numbers)." << endl;
    cout << "                                  If 0 or no extra argument is given a branch" << endl;
    cout << "                                  will be selected randomly. If 'all' is given" << endl;
    cout << "                                  NNI will be performed for all branches, e.g." << endl;
    cout << "                                  --nni 4, or --nni all." << std::endl;
    cout << "--no_branch_length / -0           do not print branch lengths. If there are no" << endl;
    cout << "                                  branch lengths in input tree the default is to" << endl;
    cout << "                                  print zero length branches in the out tree." << endl;
    cout << "                                  This argument override this and print no" << endl;
    cout << "                                  branch lengths." << endl;
    cout << "--null_short_branches [value]     set branches with shorter than given value to 0" << endl;
    cout << "--number_of_taxa / -n             get the number of tips/taxa in the tree." << endl;
    cout << "--outgroup_root / -o [taxa]	    root using most recent common ancestor of given" << endl;
    cout << "                                  taxa, e.g. -o taxa1,taxa2." << endl;
    cout << "--output [newick/nexus]           give tree format for output, nexus (nex or x" << endl;
    cout << "                                  for short), newick (new or w for short), or" << endl;
    cout << "                                  svg e.g. --output x. (default w). For svg" << endl; 
    cout << "                                  extra graphical commands can be given after a" << endl;
    cout << "                                  colon (:). Each command should be separated by" << endl;
    cout << "                                  &, and commands and arguments should be" << endl;
    cout << "                                  separated by colon. Possible commands are:" << endl;
    cout << "                                  'width' set width of figure; 'height' set" << endl;
    cout << "                                  hight of figure; 'offset' set offset between" << endl;
    cout << "                                  tips and tip label; 'strokewidth' set the" << endl;
    cout << "                                  width of the branches; 'fontsize' sets the" << endl;
    cout << "                                  size of the font used; 'font' sets which font" << endl;
    cout << "                                  to use; and 'tipcolor' sets the color of the" << endl;
    cout << "                                  tip labels given in parenthesis directly" << endl;
    cout << "                                  behind the color. 'width' and 'height' are" << endl;
    cout << "                                  mandatory if you want to set any other" << endl;
    cout << "                                  parameter than tip color. E.g. --output 'svg:" << endl;
    cout << "                                  width:300&height:400&tipcolor:red(taxon1," << endl;
    cout << "                                  taxon2,taxon3)green(taxon4)'." << endl;
    cout << "--patristic_distances / -p [sep.] get the total patristic distance to all other" << endl;
    cout << "                                  taxa in the tree for each taxon, the separator" << endl;
    cout << "                                  between different taxa, and the separator" << endl;
    cout << "                                  between taxon name and value can be specified" << endl;
    cout << "                                  (separated by colon) e.g. -p \",: | \"" << endl;
    cout << "                                  (default is new line and space). A separator" << endl;
    cout << "                                  between output from different trees can be" << endl;
    cout << "                                  given after another colon." << endl;
    /// Secret function ///
/*    cout << "--rand_node_depths                assign random node depths according to the Yule" << endl;
    cout << "                                  or Coalescence model (given as extra argument," << endl;
    cout << "                                  default Yule). The Yule model will have rate 1" << endl;
    cout << "                                  and the Coalescence branch lengths will be" << endl;
    cout << "                                  scaled by population size and generation time." << endl;
    cout << "                                  e.g. --rand_node_depths Coalescence." << endl; */
    cout << "--random_tree / -r                get a random topology (no branch lengths) with" << endl;
    cout << "                                  given number of taxa, e.g. -r 20 (default 0)." << endl;
    cout << "                                  Number of random trees can be given behind a" << endl;
    cout << "                                  colon (:), e.g. -r 20:100." << endl;
    cout << "--read_figtree_annotations        will read annotations in FigTree/treeanotator" << endl;
    cout << "                                  format (e.g. [&rate=1.0,height=3.0])." << endl;
    cout << "--relaxed_outgroup_root [taxa]    will root on the group defined as for" << endl;
    cout << "                                  --get_relaxed_outgroup." << endl;
    cout << "--set_branch_lengths / -b [value] set all branches in the tree to the given" << endl;
    cout << "                                  value, e.g. 0.5 (default 1.0)." << endl;
////////////////////
    cout << "--split [method:criterion]        splits tree based on the longest branch" << endl;
    cout << "                                  (longest_branch/l) or the mid point" << endl;
    cout << "                                  (mid_point/m) until a stop criterion set by" << endl;
    cout << "                                  --split_stop is reached. Which derived tree to" << endl;
    cout << "                                  split in each iteration can be set after :." << endl;
    cout << "                                  Either the tree with the longest branch (l;" << endl;
    cout << "                                  default for longest_branch split) or the tree" << endl;
    cout << "                                  with most tips (n; default for mid_point" << endl;
    cout << "                                  split)." << endl;
    cout << "--split_stop [stop_crit.:integer] sets criterion for when to stop splitting" << endl;
    cout << "                                  trees, either at a maximum number of trees" << endl;
    cout << "                                  (max_tree_number/t) or when all trees have" << endl;
    cout << "                                  fewer than a certain number of tips" << endl;
    cout << "                                  (max_tree_size/s). The number should be given" << endl;
    cout << "                                  togather with the specific criterion after :." << endl;
/////////////////////////
    cout << "--rooted / -R                     sets if the tree should be considered as" << endl;
    cout << "                                  rooted or not (only matters when splitting" << endl;
    cout << "                                  trees)." << endl;
    cout << "--sum_branch_length / -s          get the sum of the branch lengths in the tree" << endl;
    cout << "                                  (including root branch if length for this is" << endl;
    cout << "                                  given)." << endl;
    cout << "--verbose / -v                    get additional output." << endl;
    cout << endl;
}

