/********************************************************************
Copyright (C) 2015 Martin Ryberg

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
#ifndef FILE_PARSER
#define FILE_PARSER

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <bitset>
#include "constants.h"
#include "character_vector.h"

using namespace std;

class file_parser {
    private:
    //istream infile;
    char file_type;
    //void set_file_type( );
    void pars_relaxed_phylip( vector<character_vector>& matrix, map<char,bitset<SIZE> >& alphabet);
    bool is_whitespace(const char c) { return (c==' ' || c=='\t' || c=='\n' || c=='\r'); }
    public:
    istream* file_stream;
    file_parser (istream* input): file_type('0'), file_stream(input){};
    file_parser (): file_type('0'), file_stream(0){};
    void set_input_stream( istream* stream) { file_stream = stream; };
    //void set_data_stream( istream& stream) { datafile = stream; };
    bool set_file_type (const char* type);
    bool set_file_type();
    bool move_to_next_X_block(const char X);
    char move_to_next_block();
    char read_next_nexus_command();
    void read_translate_parameters( map<string,string>& translations );
    bool read_mrca( map<string,string>& taxonset );
    bool read_fixage( map<string,string>& ages);
    void move_to_end_of_command() {
	char character('0');
	while (*file_stream && character!=';') file_stream->get();
    }
    bool move_to_start_of_tree ();
    bool test_file_type ( const char* type ) {
	if (strcmp(type,"nexus")==0 && file_type=='N') return true;
	else if (strcmp(type,"newick")==0 && file_type=='n') return true;
       	else if (strcmp(type,"phylip")==0 && file_type=='p') return true;
	else if (strcmp(type,"fasta")==0 && file_type=='f') return true;
	else return false;
    };
    const char* get_file_type() {
	#ifdef DEBUG
	cerr << "Type as in object " << file_type << endl;
	#endif //DEBUG
	if (file_type=='N') return "nexus";
	else if (file_type=='n') return "newick";
	else if (file_type=='p') return "phylip";
	else if (file_type=='f') return "fasta";
	else return "non";
    };
    void pars_data ( vector<character_vector>& matrix, map<char,bitset<SIZE> >& alphabet) {
	if (file_type == 'p') pars_relaxed_phylip(matrix,alphabet);
    };
    bool move_to_begining_of_file() {
	if (file_stream == &cin) return false;
	else {
	    if (file_stream->eof()) file_stream->clear();
	    if (file_stream->good()) {
		file_stream->seekg(0,ios::beg);
		return true;
	    }
	}
	return false;
    };
};

#endif //FILE_PARSER
