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

#ifndef CHARACTER_VECTOR
#define CHARACTER_VECTOR

#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <vector>
#include "constants.h"
using namespace std;

class parsimony_character_vector;

class character_vector {
public:
    void add_character (bitset<SIZE> character) {
	#ifdef DEBUG
	cerr << "Trying to add " << character << endl;
	#endif // DEBUG
	characters.push_back(character);
    };
    void add_character (const unsigned int pos, bitset<SIZE>& character) {
	#ifdef DEBUG
	cerr << "Trying to add " << character << " at pos " << pos << endl;
	#endif // DEBUG
	unsigned int vector_size = characters.size();
	if (pos == vector_size) characters.push_back(character);
	else if (pos < vector_size) characters[pos] = character;
	else {
	    while (pos > vector_size) {
		characters.push_back(bitset<SIZE>());
		++vector_size;
	    }
	    characters.push_back(character);
	}
    };
    bitset<SIZE>& get_character (unsigned int pos) {
	unsigned int vector_size = characters.size();
	if (pos < vector_size) return characters[pos];
	else {
	    while (pos >= vector_size) {
		characters.push_back(bitset<SIZE>());
		++vector_size;
	    }
	    return characters.back();
	}
    };
    void reset_char() { characters.clear(); };
    void reset_taxon() { taxon.clear(); };
    void reset() { characters.clear(); taxon.clear(); };
    void clear () { reset(); };
    bool empty() { return characters.empty() && taxon.empty(); };
    unsigned int n_char() { return characters.size(); };
    unsigned int max_n_char() { return characters.max_size(); };
    string& get_taxon() { return taxon; };
    void set_taxon ( const string& name ) { taxon = name; };
    unsigned int highest_char_state() { return highest_char_state(0,characters.size()-1); }
    unsigned int highest_char_state(unsigned int start, unsigned int end) {
	unsigned int higest_state(0);
	for (unsigned int i=start; i <= end && i < characters.size(); ++i)
	    for (unsigned int j=higest_state; j < SIZE; ++j)
		if (characters[i].test(j)) higest_state = j;
	return higest_state;
    };
    void set_empty_to( bitset<SIZE> value ) {
	for (vector<bitset<SIZE> >::iterator i=characters.begin(); i != characters.end(); ++i)
	    if (i->none()) *i |= value;
    }
protected:
    vector<bitset<SIZE> > characters;
    string taxon;
    friend class parsimony_character_vector;
};

class parsimony_character_vector : public character_vector {
    public:
	parsimony_character_vector() : score(0) {};
	parsimony_character_vector(character_vector other) : score(0) {
	    characters = other.characters;
	    taxon = other.taxon;
	}
	void set_score(const unsigned int s) {
	    score = s;
	};
	unsigned int get_score () { return score; };
    private:
	unsigned int score;
};
#endif //CHARACTER_VECTOR
