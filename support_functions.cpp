#ifndef SUPPORT_FUNCTIONS
#define SUPPORT_FUNCTIONS

#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <map>
#include "constants.h"

using namespace std;

namespace support_functions {
    bitset<SIZE> pick_a_random_true_bit(const bitset<SIZE>& trait) {
	bitset<SIZE> return_bit;
	unsigned int n = trait.count();
	if (n>1) {
	    unsigned int pick = rand() % n;
	    for (unsigned int i=0; i<SIZE; ++i) {
		if (trait[i]) {
		    if (!pick) {
			return_bit.set(i);
			break;
		    }
		    --pick;
		}
	    }
	    return return_bit;
	}
	else return trait;
    }
    string pars_ARGV_string ( char [] input ) {
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
		while (fileinput) {
		    char in;
		    fileinput >> in;
		    argv_string += in;
		}
	    }
	}
	return argv_string;
    }
}
#endif //SUPPORT_FUNCTIONS
