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
#include "clustertree.h"

using namespace std;

void help ();

int main (int argc, char *argv []) {
    char method = '0'; //what method to use for clustering, default branch length over certain distance
    float cut_off = 1.00;
//    unsigned int max_size = 1;
    unsigned int name_position = 1;
    char separator = ' ';
    string file_name;
    ifstream infile;
    istream* instream = &std::cin;
    #ifdef DATABASE
    string database_data;
    #endif /*DATABASE*/
    bool quiet = false;
    for (int i=1; i < argc; ++i) {
        if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--long_branch")) method = 'l';
        else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--tip_name")) {
            method = 't';
            if (argv[i+1][0] != '-') {
                separator = argv[++i][0];
                if (argv[i+1][0] != '-') {
                    name_position = atoi(argv[++i]);
                }
            }
        }
        #ifdef DATABASE
        else if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--database")) {
            method = 'd';
            if (argv[i+1][0] != '-') database_data = argv[++i];
        }
        #endif /*DATABASE*/
        else if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"--branch_length")) method = 'b';
        else if (!strcmp(argv[i],"-c") || !strcmp(argv[i],"--cut_off")) cut_off = atof(argv[++i]);
//        else if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--max_size")) max_size = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-q") || !strcmp(argv[i],"--quiet")) quiet = true;
	else if (!strcmp(argv[i],"-f") || !strcmp(argv[i],"--file")) {
	    if (i < argc-1 && argv[i+1][0] != '-') file_name = argv[++i];
	    else {
		cerr << "-f / --file needs to be followed by a file name, e.g. -f file.tree." << endl;
		return 1;
	    }
	}
        else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) { help(); return 0; }
	else if (i == argc -1 && argv[i][0] != '-' && file_name.empty()) file_name = argv[++i];
	else {
	    cerr << "Do not recognize argument: " << argv[i] << endl;
	    return 1;
	}
    }
    if (method == '0') {
        if (!quiet) std::cerr << "Must give a method. The available methods are given by -h or --help" << endl;
        return 1;
    }
    #ifdef DATABASE
    string database;
    string table;
    string key;
    string column;
    if (method == 'd' && database_data.empty()) {
        if (!quiet) cerr << "If clustering using database entries the database, table, column for tree tip annotation, and column used for clustering must be given" << endl;
        if (!quiet) cerr << "    as a comma separated string, e.g. database,table,accno_column,species_column" << endl;
        return 1;
    }
    else if (method == 'd') {
        unsigned int k=0;
        //unsigned int j=0;
        for (unsigned int i=0; i < database_data.length(); ++i) {
            if (database_data[i] == ',') ++k;
            else if (k == 0) database += database_data[i];
            else if (k == 1) table += database_data[i];
            else if (k == 2) key += database_data[i];
            else if (k == 3) column += database_data[i];
        }
    }
    #endif /*DATABASE*/
    if ( method != 'l' && method != 't' && method != 'd' && method != 'b') {
        if (!quiet) cerr << argv[0] << " only accept the method b, d, l and t at the moment." << endl;
        return 1;
    }
    if (cut_off < 0) {
        if (!quiet) cerr << argv[0] << " only accept positive cut off values." << endl;
        return 1;
    }
    if (!file_name.empty()) {
	infile.open(file_name.c_str(), std::ifstream::in);
	if (infile.good())
	    instream = &infile;
	else {
	    cerr << "Could not open file: " << file_name << endl;
	    return 1;
	}
    }
    unsigned int n_read(0);
    while (1) {
	clustertree tree;
	tree.tree_file_parser( *instream ); //std::cin );
	if (tree.empty()) break;
	++n_read;
	cout << "### tree " << n_read << " ###" << endl;
	if (method == 'l') tree.br_length_clust_max_cout( cut_off );
	else if (method == 'b') tree.short_br_clust ( cut_off );
	else if (method == 't') tree.name_clust_cout( name_position, separator);
	#ifdef DATABASE
	else if (method == 'd') {
	    cout << "[" << database << ", " << table << ", " << key << ", " << column << "]" << endl;
	    tree.db_clust_cout( database.c_str(), table, key, column);
	}
	#endif /*DATABASE*/
    }
    return 0;
}
    
void help () {
    std::cout << "You called the help function. These are your options:" << endl << endl;
    std::cout << "--branch_length / -b          separate clusters by single link clustering based on phylogenetic distance." << endl;
    std::cout << "--cut_off / -c                the cut off to use when clustering." << endl;
    #ifdef DATABASE
    std::cout << "--database / -d               cluster based on annotations available in SQLite database. Need to be followed by a comma separated string" << endl;
    std::cout << "                                   with the database, table, column for tree tip annotation, and column used for clustering given, e.g." << endl;
    std::cout << "                                   -d database,table,accno_column,species_column" << endl;
    #endif /*DATABASE*/
    std::cout << "--file / -f [file_name]       give name of file, e.g. -f file.tree." << endl;
    std::cout << "--help / -h                   print this help." << endl;
    std::cout << "--long_branch / -l            returns taxa in clades on branches longer than cut off." << endl;
//    std::cout << "--max_size / -m               the maximum cluster size (the number of taxa in the tree is the default) not yet in use." << endl;
    std::cout << "--tip_name / -t               cluster taxa based on taxon annotation. Should be followed by a single character that separates different" << endl;
    std::cout << "                                   parts of the taxon name (default ' ') and a number giving which position in the name should be used" << endl;
    std::cout << "                                   for clustering, (default 1), e.g. -t '|' 5" << endl;
    std::cout << "--quiet / -q                  suppresses some error messages and output." << endl;
}
