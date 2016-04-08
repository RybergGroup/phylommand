#include "string_tree.h"

using namespace std;

void string_tree::delete_tree (node* position) {
    if (position != 0) {
	delete_tree(position->larger);
	delete_tree(position->smaller);
	delete position;
    }
}

string* string_tree::add_string (string value, node* position) {
    if (position->larger == 0 && position->smaller == 0 && position->lable.empty()) {
	position->lable = value;
	if (ordered) 
	    position->value = ++heighest_value;
	return &position->lable;
    }
    int comp_value = position->lable.compare(value);
    if (comp_value == 0) return &position->lable;
    else if (comp_value > 0) {
	if (position->larger == 0) position->larger = new_node();
	return add_string (value, position->larger);
    }
    else {
	if (position->smaller == 0) position->smaller = new_node();
	return add_string (value, position->smaller);
    }
}

string_tree::node* string_tree::find_string( node* position, const string query ) {
    if (position == 0) return 0;
    int comp_value = position->lable.compare(query);
    if (comp_value == 0) return position;
    else if (comp_value > 0) return find_string( position->larger, query );
    else return find_string( position->smaller, query );
}

bool string_tree::drop_string ( string object ) {
    if (!object.empty()) {
        node* present = start;
	node* previous = 0;
	node* move_node = 0;
	int comp_value = present->lable.compare(object);
	int prev_comp_value = 0;
	while (comp_value) {
	    previous = present;
	    if (comp_value > 0) present = present->larger;
	    else present = present->smaller;
	    prev_comp_value = comp_value;
	    if (present !=0) comp_value = present->lable.compare(object);
	    else comp_value=0;
	}
	if (present && previous) {
	    node* new_daughter = 0;
	    if ( present->larger !=0 ) {
		new_daughter = present->larger;
		if (present->smaller != 0) {
		    move_node = present->smaller;
		}
	    }
	    else if ( present->smaller != 0 ) {
		new_daughter = present->smaller;
	    }
	    if (prev_comp_value > 0) previous->larger = new_daughter;
	    else previous->smaller = new_daughter;
	}
	else if (present == start) {
	    if (present->larger == 0 && present->smaller == 0) {
		present->lable.clear();
		return true;
	    }
	    else {
		if ( present->larger == 0 ) {
		    start = present->smaller;
		}
		else if ( present->smaller == 0 ) {
		    start = present->larger;
		}
		else {
		    start = present->larger;
		    move_node = present->smaller;
		}
	    }
	}
	if (move_node) {
	    node* new_path = present->larger;
	    node* new_mother = present->larger;
	    int new_comp;
	    while (new_path != 0) {
		new_comp = new_path->lable.compare(move_node->lable);
		new_mother = new_path;
		if (new_comp > 0) new_path = new_path->larger;
     		else new_path = new_path->smaller;
	    }
	    if (new_comp > 0) new_mother->larger = move_node;
	    else new_mother->smaller = move_node;
	}
	if (present) {
	    delete present;
	    return true;
	}
	else return false;
    }
    else return false;
}
