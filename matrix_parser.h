#ifndef MATRIX_PARSER
#define MATRIX_PARSER

#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <map>
#include "constants.h"
#include "character_vector.h"

using namespace std;

class matrix_parser {
public:
    matrix_parser(istream& input_file, vector<character_vector>& data_matrix, map<char,bitset<SIZE> >& code) : file(input_file), matrix(data_matrix), alphabet(code), matrix_type("relaxed_phylip") {};
    void pars () {
        if (!matrix_type.compare("relaxed_phylip")) pars_relaxed_phylip();
    };
private:
    istream& file;
    map<char,bitset<SIZE> >& alphabet;
    vector<character_vector>& matrix;
    string matrix_type;
    void pars_relaxed_phylip();
};

class alphabet_parser {
public:
    alphabet_parser(istream& input_file, map<char, bitset<SIZE> >& alphabet_map) : file(input_file), alphabet(alphabet_map), file_type("white space") {};	
    void pars () {
	if (!file_type.compare("white space")) pars_whitespace();
    }
private:
    istream& file;
    string file_type;
    map<char, bitset<SIZE> >& alphabet;
    void pars_whitespace();
};

namespace alphabet {
    void set_alphabet_binary (map<char,bitset<SIZE> >& alphapet);
    char translate_bitset (const bitset<SIZE> character, map<char,bitset<SIZE> >& alphapet); 
}
#endif //MATRIX_PARSER
