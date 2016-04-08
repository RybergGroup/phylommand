#include <stdlib.h>
#include "decisiveness.h"
#include <string.h>
#include <iostream>

using namespace std;
void help ();

int main ( int argc, char *argv[]) {
    char method = '0';
    int n = 100;
    string inputstring;
    if (argc > 1) {
        for (int i=1; i < argc; ++i) {
            if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--decisiveness")) {
                method = 'd';
                if ( i < argc-1 && argv[i+1][0] != '-') inputstring = argv[++i];
                if ( i < argc-1 && argv[i+1][0] != '-') n=atoi(argv[++i]);
            }
            else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
                help();
                return 0;
            }
            else {
                std::cerr << "Unrecognized argument " << argv[i] << ". Quitting quietly." << endl;
                return 1;
            }
        }
    }
    if (method == 'd') {
        if (inputstring.empty()) {
            std::cerr << "A string of genes is needed to calculate the decisiveness (--decisiveness/-d), e.g. -d ITS,RPB2|ITS|ITS,RPB2." << endl;
            return 1;
        }
        decisiveness stat(&inputstring);
        stat.on_random_tree( n );
        cout << "The gene sampling is decisive for " << stat.get_decisiveness() << " of the trees and " << stat.get_distinguished() << " of the branches." << endl;
    }
    return 0;
}

void help () {
    std::cout << "Superstat is a command line tool to calculate statistics for supertrees/supermatices." << endl;
    std::cout << "(c) Martin Ryberg 2013." << endl << endl;
    std::cout << "Usage:" << endl << "superstata [arguments]" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--decisiveness/-d     calculates proportion of random trees for which given gene" << endl;
    std::cout << "                          sampling is decisive and the mean proportion of branches" << endl;
    std::cout << "                          that are distinguished. The genes for each taxa are given" << endl;
    std::cout << "                          as a comma (,) separated string, the genes of each taxa" << endl;
    std::cout << "                          are separated by a bar (|). The number of random trees" << endl;
    std::cout << "                          are given by a number after the genes, e.g." << endl;
    std::cout << "                          -d 'ITS,RPB2|ITS|ITS,RPB2|RPB2|RPB2|ITS' 1000." << endl;
    std::cout << "--help / -h           print this help." << endl;
    std::cout << endl;
}
