//#ifndef MATRIX_PARSER
//#define MATRIX_PARSER

#include "matrix_parser.h"

using namespace std;

void matrix_parser::pars_relaxed_phylip() {
    char character;
    string taxon;
    bitset<SIZE> trait;
    character_vector row;
    unsigned int n_char(0);
    unsigned int n_taxa(0);
    char read_mode('n');
    unsigned int pos(0);
    while (file) {
	file.get(character);
	#ifdef DEBUG
	cerr << "Read in char: " << character << " (mode: " << read_mode << ")." << endl;
	#endif //DEBUG
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
		pos = 0;
	    }
	}
	else if (read_mode == 'C') {
	    #ifdef DEBUG
	    cerr << "Setting alphabet for pos: " << pos << endl;
	    #endif // DEBUG
	    map<char, bitset<SIZE> > alphabet = regions.get_partition_alphabet(pos);
	    #ifdef DEBUG
	    cout << "Adding " << character << " to matrix" << endl;
	    #endif //DEBUG
	    if ( !row.empty() && (character == '\n' || character == '\r')) {
		read_mode='T';
		if (row.n_char() != n_char) std::cerr << "Matrix size missmatch: " << taxon << " differes in " << (row.n_char() - n_char) << " from given number of characters." << std::endl;
		matrix.push_back(row);
		row.reset();
	    }
	    else if (!row.get_taxon().empty() && alphabet.find(character) != alphabet.end()) {
		trait |= alphabet[character];
		row.add_character(trait);
		#ifdef DEBUG
		std::cerr << "Set " << row.get_taxon() << " to " << trait << std::endl;
		#endif //DEBUG
		trait.reset();
		++pos;
	    }
	}
    }
    if (matrix.size() != n_taxa) std::cerr << "Matrix size missmatch: there is a difference of " << (matrix.size() - n_taxa) << " between given number of taxa and number of read taxa." << std::endl;
}

void matrix_parser::pars_fasta() {
    char character;
    string taxon;
    bitset<SIZE> trait;
    character_vector row;
    char read_mode('n');
    unsigned int pos(0);
    unsigned int n_taxa(0);
    while (file) {
	file.get(character);
	if (character == '>') read_mode = 't';
	else if (read_mode == 't' && (character == '\n' || character == '\r') ) {
	    if (!row.empty()) {
		matrix.push_back(row);
	    }
	    row.reset();
	    ++n_taxa;
	    if (taxon.empty()) cerr << "Warning, taxon " << n_taxa << " lack name." << endl;
	    else row.set_taxon(taxon);
	    taxon.clear();
	    pos = 0;
	    read_mode = 's';
	    #ifdef DEBUG
	    cerr << "Read " << n_taxa << " taxa." << endl;
	    #endif //DEBUG
	}
	else if (character != '\n' && character != '\r' && character != '\t') {
	    if (read_mode == 't') taxon += character;
	    else if (read_mode == 's') {
		map<char, bitset<SIZE> > alphabet = regions.get_partition_alphabet(pos);
		trait |= alphabet[character];
		row.add_character(trait);
		trait.reset();
		++pos;
	    }
	}
    }
    if (!row.empty()) {
	matrix.push_back(row);
    }
}

void alphabet_parser::pars_whitespace() {
    #if DEBUG
    cerr << "Parsing alphabet" << endl;
    #endif //DEBUG
    char character;
    char type(' ');
    string value;
    bitset<SIZE> trait;
    char read_mode('b');
    while (file) {
	file.get(character);
	if (read_mode == 'b' && (character == ' ' || character == '\n' || character == '\r' || character == '\t')) continue;
	else if (read_mode == 'b') read_mode = 't';
	if (read_mode == 't') {
	    if ((character == ' ' || character == '\t') && type != ' ' ) read_mode = 'c';
	    else {
		type = character;
		read_mode = 'c';
		#if DEBUG
		cerr << "Parsing char: " << type << endl;
		#endif //DEBUG
	    }
	}
	else if (read_mode == 'c') {
	    if ((character == ' ' || character == '\t' || character == '\n' || character == '\r') && type != ' ' && !value.empty()) {
		unsigned int number = atoi(value.c_str());
		if (number >= SIZE) cerr << number << " is out of bound for '" << type << "', reading alphabet." << endl;
		else {
		    trait.set(number);
		}
		value.clear();
	    }
	    if ((character == '\n' || character == '\r') && type != ' ') {
		alphabet[type] = trait;
		trait.reset();
		value.clear();
		read_mode = 'b';
	    }
	    else if (character == '0' || character == '1' || character == '2' ||
		    character == '3' || character == '4' || character == '5' || character == '6' ||
		    character == '7' || character == '8' || character == '9')
			value += character;
	    else if (character == '\n' || character == '\r') read_mode = 'b';
	}
    }
}

namespace alphabet {
    void set_alphabet_binary (map<char,bitset<SIZE> >& alphapet) {
	if (SIZE > 1) {
	    alphapet['0'][0] = 1;
	    alphapet['1'][1] = 1;
	    alphapet['-'][0] = alphapet['-'][1] = 1;
	}
    }

    char translate_bitset (const bitset<SIZE> character, map<char,bitset<SIZE> >& alphapet) {
	for ( map<char,bitset<SIZE> >::const_iterator i = alphapet.begin(); i != alphapet.end(); ++i ) {
	    if (i->second == character) return i->first;
	}
	return '!';
    }
}
//#endif //MATRIX_PARSER
