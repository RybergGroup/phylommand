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
#include "tree.h"

using namespace std;

void help ();
void longest_branch( tree tree, bool rooted, int max_size );

const int max_trees = 1000;

main (int argc, char *argv []) {
    char method = '0';
    char split_criteria = 'n';
    char split_stop = 's';
    int max_size = 0;
    bool rooted = false;
    bool int_lables = true; // if lables on internal nodes should be printed
    tree* trees[max_trees];
    for (int i=0; i<max_trees;++i) trees[i]=0;
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-l") || !strcmp(argv[i],"--longest_branch")) {
                method = 'l';
                if ( i < argc-1 && argv[i+1][0] != '-') split_criteria = argv[++i][0];
                else split_criteria = 'l';
            }
            else if (!strcmp(argv[i],"-m") || !strcmp(argv[i],"--mid_point")) {
                method = 'm';
                if ( i < argc-1 && argv[i+1][0] != '-') split_criteria = argv[++i][0];
                else split_criteria = 'n';
            }
            else if (!strcmp(argv[i],"-t") || !strcmp(argv[i],"--max_tree_number")) {
                split_stop = 't';
                if ( i < argc-1 && argv[i+1][0] != '-') max_size = atoi(argv[++i]);
            }
            else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--max_tree_size")) {
                split_stop = 's';
                if ( i < argc-1 && argv[i+1][0] != '-') max_size = atoi(argv[++i]);
            }
            else if (!strcmp(argv[i],"-i") || !strcmp(argv[i],"--print_int_label")) {
                if ( i < argc-1 && (argv[i+1][0] == 'n' || argv[i+1][0] == 'N')) {
                    int_lables = false;
                    i++;
                }
                else if ( i < argc-1 && (argv[i+1][0] == 'y' || argv[i+1][0] == 'Y')) {
                    int_lables = true;
                    i++;
                }
                else std::cerr << "-i/--print_int_label should be followed by y(es) or n(o). Keeping default to print internal node labels." << endl;
            }
            else if (!strcmp(argv[i],"-r") || !strcmp(argv[i],"--rooted")) rooted=true;
            else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) help();
            else {
                std::cerr << "Unrecognized argument: " << argv[i] << ". Quitting." << endl;
                return 1;
            }
        }
    }
    trees[0]=new tree;
    (*trees[0]).tree_file_parser( std::cin );
    if (method == 'l' || method == 'm') {
        int last=0;
        int split=0;
        bool cont = true;
        if (max_size==0) max_size = (*trees[0]).n_tips() + 100;
        while(cont && last<max_trees) {
            cont=false;
            trees[++last] = new tree;
            if (method == 'l') {
                if (rooted) /**trees[last]=*/(*trees[split]).split_tree_at_longest_branch(trees[last]);
                else /**trees[last]=*/(*trees[split]).split_tree_at_longest_branch_unrooted(trees[last]);
            }
            else if (method == 'm') (*trees[split]).split_tree_at_midpoint(trees[last]);
            int max=0;
            double max_br=0;
            for (int i=0;i<last;++i) {
                int n=(*trees[i]).n_tips();
                if (split_stop == 's')
                    if (n > max_size) cont = true;
                else if (split_stop == 't')
                    if (last <= max_size) cont = true;
                if (split_criteria == 'n') {
                    if (n > max) {
                        max=n;
                        split=i;
                    }
                }
                else if (split_criteria == 'l') {
                    double br=(*trees[i]).longest_branch();
                    if (br > max_br) {
                        max_br=br;
                        split=i;
                    }
                }
            }
        }
        for (int i=0;i<=last;++i) {
            if (int_lables) (*trees[i]).print_newick();
            else (*trees[i]).print_newick(true,false);
            //std::cout << endl << endl;
            (*trees[i]).~tree();
            //delete trees[i];
        }
    }
    else { std::cerr << "No recognized method." << endl; }
}

void help ( ) {
    std::cout << "This is treespliter, no help is currently available. Wait for future releases" << endl;
}
