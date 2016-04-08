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
}
#endif //SUPPORT_FUNCTIONS
