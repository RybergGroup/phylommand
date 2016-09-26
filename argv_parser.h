#ifndef ARGV_PARSER
#define ARGV_PARSER

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>

using namespace std;

namespace argv_parser {
    string pars_string (const char input [] );
    void pars_sub_args (const char* argument, const char separator, vector<string>& container);
   // void pars_semicolon(char* argument, vector<string> container);
};
#endif //ARGV_PARSER
