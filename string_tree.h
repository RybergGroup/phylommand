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

#include <iostream>
#include <string.h>

using namespace std;

class string_tree {
    public:
    // Structs
    struct node {
	string lable;
	unsigned int value;
	node* larger;
	node* smaller;
    };
    // Constructors
    string_tree () {
	start = new_node();
	ordered = true;
	heighest_value=0;
    };
    string_tree (string value) {
        string_tree ();
        add_string(value);
    };
    string_tree (string value, bool type) {
	string_tree ();
	ordered = type;
	add_string (value);
    };
    //Destructor
    ~string_tree () {
	delete_tree(start);
    };
    // Interface functions
    bool is_ordered () { return ordered; }; // give if string_tree is ordered
    string* add_string (string value) { // add a string
	return add_string(value, start);	
    };
    string* find_string ( const string query ) {
	node* string_node = find_string(start, query);
	return &string_node->lable;
    };
    bool drop_string ( string object );
    private:
    unsigned int heighest_value;
    bool ordered;
    void delete_tree (node* position);
    node* start;
    node* new_node () {
	node* create = new node;
	create->larger = 0;
	create->smaller = 0;
	return create;
    };
    string* add_string (string value, node* position);
    node* find_string( node* position, const string query );

};
