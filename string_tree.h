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
