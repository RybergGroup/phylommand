#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>

using namespace std;

class indexedfasta {
    public:
    bool is_open() { return file_stream->good() && !index.empty(); };
    void open( const string& name );
    float get_comp_value(const string accno) { return float(index[accno].nonN-index[accno].N); }; // return number of non N char in float format
    void get_taxon_string(const string& accno, string& taxon) {
	taxon = index[accno].taxon_string;
	#ifdef DEBUG
	cerr << "Taxon string returned: " << taxon << endl;
	#endif // DEBUG
    };
    bool initiate_sequence_retrieval() { 
	#ifdef DEBUG
	cerr << "Initiating sequence retrieval." << endl;
	#endif //DEBUG
	seq1 = index.begin();
	return seq1 != index.end(); // returns false if index is empty
    };
    bool set_seq1 ( string accno ) { seq1 = index.find(accno); if (seq1 == index.end()) return false; else return true;};
    bool seq2_is_last();
    bool next_seq1() {
	++seq1;
	if (seq1 != index.end()) return true;
	else return false;
    };
    bool next_seq2() {
	++seq2;
	if (seq2 != index.end()) return true;
	else return false;
    };
    void set_seq2_to_seq1() { seq2 = seq1; }
    string get_accno1() { return seq1->first; };
    string get_accno2() { return seq2->first; };
    void get_sequence1 ( string& sequence ) { get_sequence (seq1, sequence); };
    void get_sequence2 ( string& sequence ) { get_sequence (seq2, sequence); };
    bool seq1_is_last_seq() {
	map<string,set>::iterator temp = seq1;
	++temp;
	if (seq1 == index.end() || temp == index.end()) return true;
	else return false;
    };
    bool seq2_is_last_seq() {
	map<string,set>::iterator temp = seq2;
	++temp;
	if (seq2 == index.end() || temp == index.end()) return true;
	else return false;
    };
    void assign_taxon_string (const string& taxon, const vector<string>& accnos);
    private:
    struct set {
	string taxon_string;
	unsigned int N;
	unsigned int nonN;
	streampos pos;
    };
    map <string, set> index;
    map<string,set>::iterator seq1;
    map<string,set>::iterator seq2;
    ifstream file;
    istream* file_stream;
    stringstream cin_holder;
    void get_sequence( map<string,set>::iterator seq, string& sequence);
};
