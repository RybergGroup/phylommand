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

#include "file_parser.h"

using namespace std;

bool file_parser::set_file_type( ) {
    if (file_stream == &cin) return false;
    #ifdef DEBUG
    cerr << "Input stream is from file" << endl;
    #endif //DEBUG
    bool return_value = false;
    streampos start = file_stream->tellg();
    #ifdef DEBUG
    cerr << "File start pos: " << file_stream->tellg() << endl;
    #endif //DEBUG
    locale loc;
    string word;
    *file_stream >> word;
    if (!word.empty()) {
	for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
	//std::cerr << word << endl;
	if (!word.compare("#nexus")) { file_type = 'N'; return_value = true; }
	else if (word[0]=='>') { file_type = 'f'; return_value = true; }
	else if (word[0]=='(') { file_type = 'n'; return_value = true; }
	else if (word[0] == '1' || word[0] == '2' || word[0] == '3' || word[0] == '4' ||
		word[0] == '5' || word[0] == '6' || word[0] == '7' || word[0] == '8' ||
		word[0] == '9') { file_type = 'p'; return_value = true; }
    }
    file_stream->seekg(start);
    #ifdef DEBUG
    cerr << "File pos when finished guessing file format: " << file_stream->tellg() << endl;
    #endif //DEBUG
    return return_value;
}

bool file_parser::set_file_type (const char* type) {
    if (strcmp(type,"nexus")==0) { file_type = 'N'; return true; }
    else if (strcmp(type,"newick")==0) { file_type = 'n'; return true; }
    else if (strcmp(type,"phylip")==0) { file_type = 'p'; return true; }
    else if (strcmp(type,"fasta")==0) { file_type = 'f'; return true; }
    else return false;
}

char file_parser::move_to_next_block(){
    if (file_type != 'N') return nexus_block::ERROR;
    locale loc;
    string word;
    while (*file_stream) {
	*file_stream >> word;
	if (!word.empty()) {
	    for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
	    if (!word.compare("begin")) {
		*file_stream >> word;
		for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
		if (!word.compare("trees;")) return nexus_block::TREES;
		else if (!word.compare("data;")) return nexus_block::DATA;
		else if (!word.compare("taxa;")) return nexus_block::TAXA;
		else if (!word.compare("phylommand;")) return nexus_block::PHYLOMMAND;
		else return nexus_block::UNKNOWN;
	    }
	}
    }
    return nexus_block::NON;
}

bool file_parser::move_to_next_X_block ( const char X ) {
    if (file_type != 'N') return false;
    streampos start = file_stream->tellg();
    char block('b');
    while (block != nexus_block::NON && block != nexus_block::ERROR) {
	#ifdef DEBUG
	cerr << "Looking for " << X << endl;
	#endif //DEBUG
	block = move_to_next_block();
	if (block == X) return true;
    }
    if (file_stream != &cin) {
	if (file_stream->eof()) file_stream->clear(ios::eofbit);
	file_stream->seekg(start);
    }
    return false;
}

char file_parser::read_next_nexus_command( ) {
    #ifdef DEBUG
    cerr << "Look for next command " << get_file_type() << endl;
    #endif //DEBUG
    if (file_type != 'N') return nexus_command::NON;
    locale loc;
    string word;
    while (*file_stream) {
	*file_stream >> word;
	#ifdef DEBUG
	cerr << "Reading word: " << word << endl;
	#endif //DEBUG
	if (!word.empty()) {
	    for (string::size_type i=0; i<word.length(); ++i) word[i]=tolower(word[i],loc);
	    if (!word.compare("translate"))
		return nexus_command::TRANSLATE;
	    else if (!word.compare("tree"))
		return nexus_command::TREE;
	    else if (!word.compare("taxset")) 
		return nexus_command::TAXSET;
	    else if (!word.compare("mrca")) 
		return nexus_command::MRCA;
	    else if (!word.compare("fixage")) 
		return nexus_command::FIXAGE;
	    else if (!word.compare("end;"))
		return nexus_command::END;
	}
    }
    return nexus_command::NON;
}

bool file_parser::read_mrca( map<string,string>& taxonset ) {
    string setname;
    string included_taxa;
    string taxon;
    char mode('n');
    char character = file_stream->get();
    while (*file_stream) {
	//if (character==' ' && character=='\t' && character=='\n' && character=='\r') {
	#ifdef DEBUG
    	cerr << character << endl;
	#endif //DEBUG
	if (is_whitespace(character) || character == ';') {
	    #ifdef DEBUG
	    cerr << "Found whitspace" << endl;
	    #endif //DEBUG
	    if (!taxon.empty()) {
		if (!included_taxa.empty()) included_taxa+=',';
		map<string,string>::const_iterator add_set = taxonset.find(taxon);
		if (add_set != taxonset.end())
		    included_taxa += add_set->second;
		else
		    included_taxa += taxon;
		taxon.clear();
	    }
	    else if (!setname.empty()) mode = 's';
	    if (character == ';') break;
	}
	else if (character == '=') mode = 's';
	else if (mode == 'n') setname += character;
	else
	    taxon += character;
	character = file_stream->get();
    }
    if (!setname.empty() && !included_taxa.empty()) {
	taxonset[setname]=included_taxa;
	return true;
    }
    else return false;
}

bool file_parser::read_fixage( map<string,string>& ages) {
    #ifdef DEBUG
    cerr << "Reading fixage." << endl;
    #endif //DEBUG
    string setname;
    string word;
    string age;
    char mode('S');
    char character = file_stream->get();
    while (*file_stream) {
	#ifdef DEBUG
	cerr << character << endl;
	#endif //DEBUG
	if (is_whitespace(character) || character == '=' || character == ';') {
	    if (!word.empty()) {
		if (mode == 'S' && !word.compare("taxon")) mode='n';
		else if (mode == 'n') {
		    setname = word;
		    mode='S';
		}
		else if (mode == 'S' && !word.compare("age")) mode = 'a';
		else if (mode == 'a') {
		    age = word;
		    mode = 'S';
		}
		#ifdef DEBUG
		cerr << word << ' ' << mode << endl;
		#endif //DEBUG
		word.clear();
	    }
	    if (character == ';') break;
	}
	else
	    word += character;
	character = file_stream->get();
    }
    if (!setname.empty() && !age.empty()) {
	ages[setname]=age;
	return true;
    }
    else return false;
}

void file_parser::read_translate_parameters( map<string,string>& translations ) {
    if (file_type != 'N') return;
    #ifdef DEBUG
    cerr << "Reading translation block." << endl;
    #endif //DEBUG
    string key;
    string value;
    char type = 'k';
    char temp;
    temp = file_stream->get();
    unsigned int number_of_sq_brakets(0);
    while (*file_stream) {
	if (temp == '[') {
	    ++number_of_sq_brakets;
	    while (*file_stream && number_of_sq_brakets) {
		#ifdef DEBUG
		cerr << "comment: " << temp << endl;
		#endif //DEBUG
		temp = file_stream->get();
		if (temp == ']') --number_of_sq_brakets;
		else if (temp == '[') ++number_of_sq_brakets;
	    }
	}
	else if (type == 'k' || type == 'K') {
	    if (temp != ' ' && temp != '\n' && temp != '\t' && temp != '\r') {
		key += temp;
		type = 'K';
	    }
	    else if (type == 'K') type = 'v';
	}
	else if (type == 'v') {
	    if (temp == ',') {
		type = 'k';
		while (!value.empty() && value.back() == ' ') value.pop_back();
		translations[key]=value;
		key.clear(); value.clear();
	    }
	    else if (temp == '\'' || temp == '"') {
		char delimiter = temp;
		temp = file_stream->get();
		while (*file_stream && temp != delimiter) {
		    value += temp;
		    temp = file_stream->get();
		}
		continue;
	    }
	    else if (!(temp == ' ' && value.empty()) && temp != '\n' && temp != '\t' && temp != '\r') value += temp;
	}
	temp = file_stream->get();
	if (temp == ';') {
	    while (!value.empty() && value.back() == ' ') value.pop_back();
	    translations[key]=value;
	    break;
	}
    }
    #ifdef DEBUG
    cerr << "Read: " << translations.size() << " taxon translations." << endl;
    #endif //DEBUG
}

bool file_parser::move_to_start_of_tree ( ) {
    if (file_type != 'N') return false;
    char temp = file_stream->get();
    while (*file_stream && temp != '=') temp = file_stream->get();
    if (temp == '=') return true;
    else return false;
}

void file_parser::pars_relaxed_phylip( vector<character_vector>& matrix, map<char,bitset<SIZE> >& alphabet) {
    char character;
    string taxon;
    bitset<SIZE> trait;
    character_vector row;
    unsigned int n_char(0);
    unsigned int n_taxa(0);
    char read_mode('n');
    while (*file_stream) {
        file_stream->get(character);
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
            }
        }
        else if (read_mode == 'C') {
            if (character == '\n') {
                read_mode='T';
                if (row.n_char() != n_char) std::cerr << "Matrix size missmatch: " << taxon << " differes in " << (row.n_char() - n_char) << " from given number of characters." << std::endl;
                matrix.push_back(row);
                row.reset();
            }
            else if (!row.get_taxon().empty() && alphabet.find(character) != alphabet.end()) {
                trait |= alphabet[character];
                row.add_character(trait);
                #ifdef DEBUG
                std::cerr << "Set " << taxon << " to " << trait << std::endl;
                #endif //DEBUG
                trait.reset();
            }
        }
    }
    if (matrix.size() != n_taxa) std::cerr << "Matrix size missmatch: there is a difference of " << (matrix.size() - n_taxa) << " between given number of taxa and number of read taxa." << std::endl;
}

