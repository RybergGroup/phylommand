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
    string pars_ARGV_string ( char [] input );
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
