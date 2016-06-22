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

class partitions {
public:
    void add_alphabet( string name, map<char, bitset<SIZE> > alphabet) { alphabets[name] = alphabet; };
    void add_partition( unsigned int start, unsigned int end, const string& name, const string& alphabet) {
	values.push_back(partition(start,end,name,get_alphabet_iterator_by_name(alphabet)));
    };
    map<char, bitset<SIZE> > get_partition_alphabet(unsigned int pos) {
	vector<partition>::iterator part = get_partition(pos);
	if (part == values.end()) {
	    get_partition_alphabet("default");
	}
	return part->get_alphabet();
    };
    map<char, bitset<SIZE> > get_partition_alphabet(const string& name) {
	vector<partition>::iterator part = get_partition(name);
	return part->get_alphabet();
    };
private:
    class partition {
    public:
	partition(unsigned int begining, unsigned int ending, const string& type, map<string, map<char, bitset<SIZE> > >::const_iterator code) : start(begining), end(ending), name(type), alphabet(code) { };
	bool in_partition(unsigned int pos) {
	    return pos >= start && pos <= end;
	};
	bool is_named(const string& nomum) { return !name.compare(nomum); };
	map<char, bitset<SIZE> > get_alphabet() { return alphabet->second; };
    private:
	unsigned int start;
	unsigned int end;
	string name;
	map<string, map<char, bitset<SIZE> > >::const_iterator alphabet;
    };
    map<string, map<char, bitset<SIZE> > >::const_iterator get_alphabet_iterator_by_name ( const string& name) { return alphabets.find(name); };
    vector<partition> values;
    map<string, map<char, bitset<SIZE> > > alphabets;
    vector<partition>::iterator get_partition ( unsigned int pos ) {
	for (vector<partition>::iterator i=values.begin(); i != values.end(); ++i) {
	    if (i->in_partition(pos)) return i;
	}
    };
    vector<partition>::iterator get_partition ( const string& name ) {
        for (vector<partition>::iterator i=values.begin(); i != values.end(); ++i) {
            if (i->is_named(name)) return i;
        }
    };

};

class matrix_parser {
public:
    matrix_parser(istream& input_file, vector<character_vector>& data_matrix, partitions& codes) : file(input_file), matrix(data_matrix), regions(codes), matrix_type("relaxed_phylip") {};
    void pars () {
        if (!matrix_type.compare("relaxed_phylip")) pars_relaxed_phylip();
    };
private:
    istream& file;
    partitions regions;
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
};
#endif //MATRIX_PARSER
