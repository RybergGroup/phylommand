#ifndef MATRIX_PARSER
#define MATRIX_PARSER

#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <vector>
#include "constants.h"

class matrix_parser {
public:
    matrix_parser(istream& input_file, vector<character_vector>& data_matrix,map<char,bitset<SIZE> >& code) : file(input_file), matrix(data_matrix), alphabet(code), matrix_type("relaxed_phylip") {};
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
/*
class alphabet_parser {
public:
    
private:
    istream& file;
    map<char, bitseat<SIZE> >& alphabet;
}
*/

void matrix_parser::pars_relaxed_phylip() {
    char character;
    string taxon;
    bitset<SIZE> trait;
    character_vector row;
    unsigned int n_char(0);
    unsigned int n_taxa(0);
    char read_mode('n');
    while (file) {
	file.get(character);
	if (read_mode == 'n') {
	    if (character == ' ' && n_taxa > 0) read_mode = 'N';
	    if (character != ' ') {
		n_taxa *= 10;
		n_taxa += character - '0';
	    }
	}
	else if (read_mode == 'N') {
	    if (character == '\n') read_mode = 'T';
	    else if (character != ' ') {
		n_char *= 10;
		n_char += character - '0';
	    }
	}
	else if (read_mode == 'T') {
	    if (character != ' ') taxon += character;
	    else if (!taxon.empty() && character == ' ') {
		row.set_taxon(taxon);
		taxon.clear();
		read_mode = 'C';
	    }
	}
	else if (read_mode == 'C') {
	    if (character == '\n') {
		read_mode='T';
		if (row.n_char() != n_char) std::cerr << "Matrix size missmatch: " << taxon << " differes in " << (row.n_char() - n_char) << " from given number of characters." << std::endl;
		matrix.push_back(row);
		row.reset();
	    }
	    else if (!row.get_taxon().empty() && alphabet.find(character) != alphabet.end()) {
		trait |= alphabet[character];
		row.add_character(trait);
		#ifdef DEBUG
		std::cerr << "Set " << taxon << " to " << trait << std::endl;
		#endif //DEBUG
		trait.reset();
	    }
	}
    }
    if (matrix.size() != n_taxa) std::cerr << "Matrix size missmatch: there is a difference of " << (matrix.size() - n_taxa) << " between given number of taxa and number of read taxa." << std::endl;
}

void set_alphabet_binary (map<char,bitset<SIZE> >& alphapet) {
    alphapet['0'][0] = 1;
    alphapet['1'][1] = 1;
    alphapet['-'][0] = alphapet['-'][1] = 1;
}
#endif //MATRIX_PARSER
