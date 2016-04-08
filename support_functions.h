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
    string get_file_type( istream& infile );
    bitset<SIZE> pick_a_random_true_bit(const bitset<SIZE>& trait);
}

namespace nexus_support {
    const char TREE = 't';
    const char TRANSLATE = 'r';
    const char END = 'e';
    const char NON = 'n';
    bool move_to_next_tree_block (istream& infile );
    char read_tree_block_command( istream& infile );
    void read_translate_parameters( istream& infile, map<string,string>& translations );
    bool move_to_start_of_tree ( istream& infile );
}

#endif //SUPPORT_FUNCTIONS
