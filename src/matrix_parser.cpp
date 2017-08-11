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
	    if (character == '\n' || character == '\r') read_mode = 'T';
	    else if (character != ' ') {
		n_char *= 10;
		n_char += character - '0';
	    }
	}
	else if (read_mode == 'T') {
	    if (character == '\n') continue;
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
    if ( matrix.size() != n_taxa) std::cerr << "Matrix size missmatch: there is a difference of " << (matrix.size() - n_taxa) << " between given number of taxa and number of read taxa." << std::endl;
}

void matrix_parser::pars_fasta() {
    #ifdef DEBUG
    cerr << "Parsing fasta." << endl;
    #endif //DEBUG
    char character;
    string taxon;
    bitset<SIZE> trait;
    character_vector row;
    char read_mode('n');
    unsigned int pos(0);
    unsigned int n_taxa(0);
    while (file) {
	file.get(character);
	#ifdef DEBUG
	cerr << "Read in char: " << character << " (mode: " << read_mode << ")." << endl;
	#endif //DEBUG
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
	else if (character != '\n' && character != '\r' && character != '\t' && character != ' ') {
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
void matrix_parser::pars_nexus() {
    #ifdef DEBUG
    cerr << "Parsing nexus DATA block." << endl;
    #endif //DEBUG
    char character;
    locale loc;
    string taxon;
    bitset<SIZE> trait;
    character_vector row;
    unsigned int n_char(0);
    unsigned int n_taxa(0);
    char read_mode('0');
    unsigned int pos(0);
    while (file) {
	file.get(character);
	#ifdef DEBUG
	cerr << "Read in char: " << character << " (mode: " << read_mode << ")." << endl;
	#endif //DEBUG
	if (read_mode == '0') {
	    if (character == ' ' || character == ';' || character == '=' || character == '\n' || character == '\r') {
	       	if (!taxon.empty()) {
		    for (unsigned int i=0; i < taxon.length(); ++i) taxon[i] = tolower(taxon[i],loc);
		    if (!taxon.compare("matrix")) {
			read_mode = 'T';
			if (n_taxa == 0) cerr << "Did not find the number of taxa in the data." << endl;
			if (n_char == 0) cerr << "Did not find the number of characters in the data." << endl;
		    }
		    else if (!taxon.compare("ntax")) read_mode = 'n';
		    else if (!taxon.compare("nchar")) read_mode = 'N';
		    taxon.clear();
		}
	    }
	    else taxon += character;
	}
	else if (read_mode == 'n') {
	    if (character == ' ' || character == '\n' || character == '\r' || character == '\t' || character == ';') { if (n_taxa > 0) read_mode = '0'; }
	    else {
		n_taxa *= 10;
		n_taxa += character - '0';
	    }
	}
	else if (read_mode == 'N') {
	    if (character == ' ' || character == '\n' || character == '\r' || character == '\t' || character == ';') { if (n_char > 0) read_mode = '0'; }
	    else {
		n_char *= 10;
		n_char += character - '0';
	    }
	}
	else if (read_mode == 'T') {
	    if (character == ';') break;
	    else if (character != ' ' && character != '\n' && character != '\r' && character != '\t') taxon += character;
	    else if (!taxon.empty() && character == ' ') {
		row.set_taxon(taxon);
		taxon.clear();
		read_mode = 'C';
		pos = 0;
	    }
	}
	else if (read_mode == 'C') {
	    if (character == ';') break;
	    #ifdef DEBUG
	    cerr << "Setting alphabet for pos: " << pos << endl;
	    #endif // DEBUG
	    map<char, bitset<SIZE> > alphabet = regions.get_partition_alphabet(pos);
	    #ifdef DEBUG
	    cout << "Adding " << character << " to matrix" << endl;
	    #endif //DEBUG
	    if ( !row.empty() && (character == '\n' || character == '\r') ) {
		read_mode='T';
		if (n_char != 0 && row.n_char() != n_char) std::cerr << "Matrix size missmatch: " << taxon << " differes in " << (row.n_char() - n_char) << " from given number of characters." << std::endl;
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
    if (n_taxa != 0 && matrix.size() != n_taxa) std::cerr << "Matrix size missmatch: there is a difference of " << (matrix.size() - n_taxa) << " between given number of taxa and number of read taxa." << std::endl;
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
		if (number >= SIZE) cerr << number << " is out of bound for '" << type << "', when reading alphabet." << endl;
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
    void set_alphabet_binary (map<char,bitset<SIZE> >& alphabet) {
	if (SIZE > 1) {
	    alphabet['0'][0] = 1;
	    alphabet['1'][1] = 1;
	    alphabet['-'][0] = alphabet['-'][1] = 1;
	}
    }
    void set_alphabet_dna (map<char,bitset<SIZE> >& alphabet) {
	if (SIZE > 3) {
	    alphabet['A'].set(0);
	    alphabet['a'].set(0);
	    alphabet['G'].set(1);
	    alphabet['g'].set(1);
	    alphabet['C'].set(2);
	    alphabet['c'].set(2);
	    alphabet['T'].set(3);
	    alphabet['t'].set(3);
	    alphabet['-'].reset();
	    alphabet['R'].set(0); alphabet['R'].set(1);
	    alphabet['r'].set(0); alphabet['r'].set(1);
	    alphabet['Y'].set(2); alphabet['Y'].set(3);
	    alphabet['y'].set(2); alphabet['y'].set(3);
	    alphabet['S'].set(1); alphabet['S'].set(2);
	    alphabet['s'].set(1); alphabet['s'].set(2);
	    alphabet['W'].set(0); alphabet['W'].set(3);
	    alphabet['w'].set(0); alphabet['w'].set(3);
	    alphabet['K'].set(1); alphabet['K'].set(3);
	    alphabet['k'].set(1); alphabet['k'].set(3);
	    alphabet['M'].set(0); alphabet['M'].set(2);
	    alphabet['m'].set(0); alphabet['m'].set(2);
	    alphabet['B'].set(1); alphabet['B'].set(2); alphabet['B'].set(3);
	    alphabet['b'].set(1); alphabet['b'].set(2); alphabet['b'].set(3);
	    alphabet['D'].set(0); alphabet['D'].set(1); alphabet['D'].set(3);
	    alphabet['d'].set(0); alphabet['d'].set(1); alphabet['d'].set(3);
	    alphabet['H'].set(0); alphabet['H'].set(2); alphabet['H'].set(3);
	    alphabet['h'].set(0); alphabet['h'].set(2); alphabet['h'].set(3);
	    alphabet['V'].set(0); alphabet['V'].set(1); alphabet['V'].set(2);
	    alphabet['v'].set(0); alphabet['v'].set(1); alphabet['v'].set(2);
	    alphabet['N'].set(0); alphabet['N'].set(1); alphabet['N'].set(2); alphabet['N'].set(3);
	    alphabet['n'].set(0); alphabet['n'].set(1); alphabet['n'].set(2); alphabet['n'].set(3);
	    alphabet['.'].reset();
	    //alphabet['.'].set(0); alphabet['.'].set(1); alphabet['.'].set(2); alphabet['.'].set(3);
	    #ifdef DEBUG
	    for (map<char, bitset<SIZE> >::const_iterator i = alphabet.begin(); i != alphabet.end(); ++i) {
		cerr << i->first << " = " << i->second.to_string() << endl;
	    }
	    #endif //DEBUG
	}
    }
    void set_alphabet_amino_acid (map<char,bitset<SIZE> >& alphabet) {
	if (SIZE > 22) {
	    alphabet['A'].set(0);
	    alphabet['a'].set(0);
	    alphabet['B'].set(1);
	    alphabet['b'].set(1);
	    alphabet['C'].set(2);
	    alphabet['c'].set(2);
	    alphabet['D'].set(3);
	    alphabet['d'].set(3);
	    alphabet['E'].set(4);
	    alphabet['e'].set(4);
	    alphabet['F'].set(5);
	    alphabet['f'].set(5);
	    alphabet['G'].set(6);
	    alphabet['g'].set(6);
	    alphabet['H'].set(7);
	    alphabet['h'].set(7);
	    alphabet['I'].set(8);
	    alphabet['i'].set(8);
	    alphabet['K'].set(9);
	    alphabet['k'].set(9);
	    alphabet['L'].set(10);
	    alphabet['l'].set(10);
	    alphabet['M'].set(11);
	    alphabet['m'].set(11);
	    alphabet['N'].set(12);
	    alphabet['n'].set(12);
	    alphabet['P'].set(13);
	    alphabet['p'].set(13);
	    alphabet['Q'].set(14);
	    alphabet['q'].set(14);
	    alphabet['R'].set(15);
	    alphabet['r'].set(15);
	    alphabet['S'].set(16);
	    alphabet['s'].set(16);
	    alphabet['T'].set(17);
	    alphabet['t'].set(17);
	    alphabet['V'].set(18);
	    alphabet['v'].set(18);
	    alphabet['W'].set(19);
	    alphabet['w'].set(19);
	    alphabet['X'].set(20);
	    alphabet['x'].set(20);
	    alphabet['Y'].set(21);
	    alphabet['y'].set(21);
	    alphabet['Z'].set(22);
	    alphabet['z'].set(22);
	    alphabet['-'].reset();
	    alphabet['.'].reset();
	}
    }

    char translate_bitset (const bitset<SIZE> character, map<char,bitset<SIZE> >& alphabet) {
	for ( map<char,bitset<SIZE> >::const_iterator i = alphabet.begin(); i != alphabet.end(); ++i ) {
	    if (i->second == character) return i->first;
	}
	return '!';
    }
}
//#endif //MATRIX_PARSER
