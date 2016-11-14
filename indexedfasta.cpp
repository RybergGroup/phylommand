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

#include "indexedfasta.h"

void indexedfasta::open ( const string& name ) {
    if (name.empty()) {
	while (cin) {
	    char temp;
	    cin >> noskipws >> temp;
	    #ifdef DEBUG
	    cerr << temp;
	    #endif //DEBUG
	    cin_holder << temp;
	}
	file_stream = &cin_holder;
    }
    else {
	file.open(name.c_str());
	file_stream = &file;
    }
    string line;
    string accno;
    unsigned int seqno(0);
    while (getline(*file_stream,line)) {
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
	    index[accno] = {taxon_string, 0, 0, file_stream->tellg()};
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
    file_stream->clear();
    file_stream->seekg(0,file_stream->beg);
};

void indexedfasta::get_sequence ( map<string,set>::iterator seq, string& sequence) {
    sequence.clear();
    if (seq != index.end()) {
	file_stream->seekg(seq->second.pos,file_stream->beg);
	string line;
	while (getline(*file_stream, line)) {
	    if (line[0] != '>') sequence += line;
	    else break;
	    if (file_stream->peek() == EOF) break;
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

void indexedfasta::assign_taxon_string (const string& taxon, const vector<string>& accnos) {
    for (vector<string>::const_iterator i = accnos.begin(); i != accnos.end(); ++i) {
	map<string,set>::iterator pos = index.find(*i);
	if (pos != index.end()) pos->second.taxon_string = taxon;
	else cerr << "Warning!!! No sequence was found for '" << *i << "'" << endl;
    }
}
