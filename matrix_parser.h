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

class partitions { // class for vector of partitions
public:
    void add_alphabet( string name, map<char, bitset<SIZE> > alphabet) { alphabets[name] = alphabet; };
    void add_partition( unsigned int start, unsigned int end, const string& name, const string& alphabet) {
	values.push_back(partition(start,end,name,get_alphabet_iterator_by_name(alphabet)));
    };
    map<char, bitset<SIZE> > get_partition_alphabet(unsigned int pos) {
	#ifdef DEBUG
	cerr << "Gettin partition for: " << pos << endl;
	#endif //DEBUG
	vector<partition>::iterator part = get_partition(pos);
	if (part == values.end()) {
	    #ifdef DEBUG
	    cerr << "No partition found." << endl;
	    #endif //DEBUG
	    return get_partition_alphabet("default");
	}
	return part->get_alphabet();
    };
    map<char, bitset<SIZE> > get_partition_alphabet(const string& name) {
	vector<partition>::iterator part = get_partition(name);
	return part->get_alphabet();
    };
private:
    class partition { //class for each partition
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
	return values.end();
    };
    vector<partition>::iterator get_partition ( const string& name ) {
        for (vector<partition>::iterator i=values.begin(); i != values.end(); ++i) {
            if (i->is_named(name)) return i;
        }
	return values.end();
    };

};

class matrix_parser {
public:
    matrix_parser(istream& input_file, vector<character_vector>& data_matrix, partitions& codes) : file(input_file), regions(codes), matrix(data_matrix), matrix_type("fasta") {};
    matrix_parser(istream& input_file, vector<character_vector>& data_matrix, partitions& codes, string type) : file(input_file), regions(codes), matrix(data_matrix), matrix_type(type) {};
    void pars () {
        if (!matrix_type.compare("fasta")) pars_fasta();
        else if (!matrix_type.compare("phylip")) pars_relaxed_phylip();
	else if (!matrix_type.compare("nexus")) pars_nexus();
    };
private:
    istream& file;
    partitions regions;
    vector<character_vector>& matrix;
    string matrix_type;
    void pars_relaxed_phylip();
    void pars_fasta();
    void pars_nexus();
};

class alphabet_parser {
public:
    alphabet_parser(istream& input_file, map<char, bitset<SIZE> >& alphabet_map) : file(input_file), file_type("white space"), alphabet(alphabet_map)  {};	
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
    void set_alphabet_binary (map<char,bitset<SIZE> >& alphabet);
    void set_alphabet_dna (map<char,bitset<SIZE> >& alphabet);
    void set_alphabet_amino_acid (map<char,bitset<SIZE> >& alphabet);
    char translate_bitset (const bitset<SIZE> character, map<char,bitset<SIZE> >& alphabet); 
};
#endif //MATRIX_PARSER
