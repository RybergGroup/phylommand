#ifndef ARGV_PARSERPP
#define ARGV_PARSERPP

#include "argv_parser.h"

string argv_parser::pars_string ( char input [] ) {
    string argv_string;
    bool file_pars(false);
    for (unsigned int i(0); input[i] != '\0'; ++i) {
	if (input[i] == ':' && !argv_string.compare("file")) {
	    file_pars = true;
	    argv_string.clear();
	}
	else argv_string += input[i];
    }
    if (file_pars) {
	ifstream fileinput;
	fileinput.open(argv_string.c_str());
	argv_string.clear();
	if (fileinput.is_open()) {
	    char in;
	    while (fileinput.get(in)) {
		//fileinput.get(in);// >> in;
		if (in != '\n' && in != '\r') argv_string += in;
	    }
	    if (fileinput.is_open()) fileinput.close();
	}
    }
    return argv_string;
}

void argv_parser::pars_sub_args (char* argument, const char separator, vector<string>& container) {
    unsigned int j(0);
    bool escape(false);
    container.push_back(string());
    while (argument[j] != '\0') {
	if (argument[j] == '\\' && !escape) escape = true;
	else {
	    if (escape) {
		if (argument[j] == 'n') argument[j] = '\n';
		else if (argument[j] == 'r') argument[j] = '\r';
		else if (argument[j] == 't') argument[j] = '\t';
	    }
	    if (argument[j] == separator && !escape) container.push_back(string());
	    else container.back() += argument[j];
	    escape = false;
	}
	++j;
    }
}

#endif //ARGV_PARSER
#define ARGV_PARSER

