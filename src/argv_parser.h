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
    //void pars_clade_rates (const char* argument, rate_model& rate_mod, vector<string>& taxon_string);
    //void pars_rates_and_times(const char* argument, rate_model& rate_mod);
   // void pars_semicolon(char* argument, vector<string> container);
};
#endif //ARGV_PARSER
