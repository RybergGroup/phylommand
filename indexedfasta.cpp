#include "indexedfasta.h"

void indexedfasta::open ( const char* name ) {
    file.open(name);
    string line;
    string accno;
    unsigned int seqno(0);
    while (getline(file,line)) {
	if (line[0] == '>') {
	    ++seqno;
	    int length = line.length();
	    int mode = 'a';
	    accno.clear();
	    string taxon_string;
	    for (int i=1; i < length; ++i) {
		if (line[i] == '|' && (mode == 't' || mode == 'T')) break;
		else if (line[i] == '|' && mode == 'a') mode = 't';
		else if (mode == 'a' && line[i] != ' ') accno += line[i];
		else if (mode == 't' && line[i] != ' ') { taxon_string += line[i]; mode = 'T'; }
		else if (mode == 'T') taxon_string += line[i];
	    }
	    if (accno.empty()) {
		stringstream converter;
		converter << seqno;
		accno=converter.str();
	    }
	    index[accno] = {taxon_string, 0, 0, file.tellg()};
	    #ifdef DEBUG
	    cerr << accno << ' ' << index[accno].taxon_string << ' ' << index[accno].pos << endl;
	    #endif //DEBUG
	}
	else {
	    int length = line.length();
	    for (int i =0; i < length; ++i) {
		if (line[i] != ' ' && line[i] != '\n'  && line[i] != '\r'  && line[i] != '\t' ) {
		    if (line[i] == 'N' || line[i] == 'n') ++index[accno].N;
		    else ++index[accno].nonN;
		}
	    }
	}
    }
    file.clear();
    file.seekg(0,file.beg);
};

void indexedfasta::get_sequence ( map<string,set>::iterator seq, string& sequence) {
    sequence.clear();
    if (seq != index.end()) {
	file.seekg(seq->second.pos,file.beg);
	string line;
	while (getline(file, line)) {
	    if (line[0] != '>') sequence += line;
	    else break;
	    if (file.peek() == EOF) break;
	}
    }
}

bool indexedfasta::seq2_is_last() {
    if (index.empty()) return true;
    map<string,set>::iterator temp = index.end();
    --temp;
    if (seq2 == temp) return true;
    else return false;
}
