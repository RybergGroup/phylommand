#include "support_functions.h"

using namespace std;

string support_functions::get_file_type( istream& infile ) {
    streampos start = infile.tellg();
    locale loc;
    string word;
    infile >> word;
    string return_string = "non";
    if (!word.empty()) {
	for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
	//std::cerr << word << endl;
	if (!word.compare("#nexus")) return_string = "nexus";
	else if (word[0]=='>') return_string = "fasta";
    }
    infile.seekg(start);
    return return_string;
}
bitset<SIZE> support_functions::pick_a_random_true_bit(const bitset<SIZE>& trait) {
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

bool nexus_support::move_to_next_tree_block (istream& infile ) {
    streampos start = infile.tellg();
    locale loc;
    string word;
    while (infile) {
	infile >> word;
	if (!word.empty()) {
	    for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
	    if (!word.compare("begin")) {
		infile >> word;
		//std::cerr << word << endl;
		for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
		if (!word.compare("trees;")) return true;
	    }
	}
    }
    if (infile.eof()) infile.clear(ios::eofbit);
    if (infile != cin) infile.seekg(start);
    return false;
}
char nexus_support::read_tree_block_command( istream& infile ) {
    locale loc;
    string word;
    while (infile) {
	infile >> word;
	if (!word.empty()) {
	    for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
	    if (!word.compare("translate")) {
		return 'r';
	    }
	    if (!word.compare("tree")) {
		return 't';
	    }
	    if (!word.compare("end;")) {
		return 'e';
	    }
	}
    }
    return 'n';
}
void nexus_support::read_translate_parameters( istream& infile, map<string,string>& translations ) {
    string key;
    string value;
    char type = 'k';
    char temp;
    temp = infile.get();
    while (infile) {
	if (type == 'k' || type == 'K') {
	    if (temp != ' ' && temp != '\n' && temp != '\t' && temp != '\r') {
		key += temp;
		type = 'K';
	    }
	    else if (type == 'K') type = 'v';
	}
	else if (type == 'v') {
	    if (temp == ',') {
		type = 'k';
		translations[key]=value;
		key.clear(); value.clear();
	    }
	    else if (temp != ' ' && temp != '\n' && temp != '\t' && temp != '\r') value += temp;
	}
	temp = infile.get();
	if (temp == ';') {
	    translations[key]=value;
	    break;
	}
    }
}
bool nexus_support::move_to_start_of_tree ( istream& infile ) {
    char temp=infile.get();
    while (infile && temp != '=') temp = infile.get();
    if (temp == '=') return true;
    else return false;
}
