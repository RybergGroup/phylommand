/********************************************************************
Copyright (C) 2020 Martin Ryberg

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

#include "tree.h"

/*** Support functions for the tree class ***/
namespace support_functions {
    bitset<SIZE> pick_a_random_true_bit(const bitset<SIZE>& trait) {
        bitset<SIZE> return_bit;
        unsigned int n = trait.count();
        if (n>1) {
            unsigned int pick = rand() % n;
            for (unsigned int i=0; i<SIZE; ++i) {
                if (trait[i]) {
                    if (!pick) {
                        return_bit.set(i);
                        break;
                    }
                    --pick;
                }
            }
            return return_bit;
        }
        else return trait;
    }
}

/***********************************
*** Functions for the tree class ***
***********************************/

/** Define tree for node labels ***/
string_tree tree::nodelabels;

bool tree::random_seeded(false);
/** Function to delete tree ***
*** beginning with leaf      **/
void tree::destroy_tree (node *leaf) {
    if (leaf != 0 ) {
        destroy_tree(leaf->left); //Recursively destroy left nodes
        destroy_tree(leaf->right); //Recursively destroy right nodes
	/*if (leaf->other !=0) {
	    delete [] leaf->other;
	    leaf->other = 0;
	}*/
        delete leaf; //Actually destroy present node
    }
}

/** Function to parse tree from newik formated file                   ***
*** Basically follow a algorithm given by Revell at:                  ***
*** http://phytools.blogspot.com/2011/02/building-tree-in-memory.html **/
void tree::tree_file_parser( istream& infile, map<string,string> label_translation, bool read_extra_node_annotation ) {
    #ifdef DEBUG
    cerr << "Start pos: " << infile.tellg() << endl;
    #endif //DEBUG
    node *present_node = root; //the tree will be parsed from the root
    node *parent_node; //keep track of which node we came from
    string label_string;
    string extra_label_string;
    string branch_length_annotation;
    char read_mode = 's';
    char temp; //store input one character at the time
    infile >> temp; //read input one character at the time
    while (infile) {
//	#ifdef DEBUG
//	cerr << temp << endl;
//	#endif //DEBUG
        //ignore things in square brackets
        if ( temp == '[' ) {
	    infile >> temp;
	    if (infile.bad()) break;
	    if (temp != ']') {
		if (read_extra_node_annotation && temp == '&') { // square brachet and & means extra info
		    extra_label_string+= '[';
		    while (infile && temp != ']') {
			extra_label_string+=temp; // save it to label
			infile >> temp;
		    }
		    if (temp == ']') extra_label_string+=temp;
		}
		else {
		    int n_right_square_brachets = 1; //keep track of number of right brackets
		    while (infile && n_right_square_brachets > 0) { //while more right than left brackets have been read
			if ( temp == '[' ) ++n_right_square_brachets; //count right brackets
			if ( temp == ']' ) --n_right_square_brachets; //a left bracket neutralizes a right bracket
			infile >> temp; //keep reading input one character at the time
		    }
		    continue;
		}
            }
        }
        //if left bracket create a new left node
        else if ( temp == '(' ) {
	    read_mode = 'l';
	    //if (present_node==root) tree_comment = extra_label_string;
            present_node->left = new node; //create the left node
            parent_node = present_node; //store present node
            present_node = present_node->left; //move to new node
            present_node->parent = parent_node; //set previous node as parent node
	    extra_label_string.clear();
	    label_string.clear();
        }
        //if , go back one node and create new right node
        else if ( temp == ',' ) {
	    read_mode = 'l';
	    if (!label_string.empty() || !extra_label_string.empty()) {
		if (!label_translation.empty()) {
		    map<string,string>::const_iterator translation = label_translation.find(label_string);
		    if (translation != label_translation.end()) label_string = translation->second;
		}
		present_node->nodelabel = nodelabels.add_string(label_string + extra_label_string);
		label_string.clear();
		extra_label_string.clear();
		#ifdef DEBUG
		cerr << *present_node->nodelabel << endl;
		#endif //DEBUG
	    }
	    if (!branch_length_annotation.empty()) {
		present_node->branchlength = atof(branch_length_annotation.c_str());
		branch_length_annotation.clear();
	    }
            present_node = present_node->parent; //go back one node
            if (present_node->right == 0) {  // if not already a branch there
                present_node->right = new node; //create new right node
                parent_node = present_node; //store present node location
                present_node = present_node->right; //go to newly created node
                present_node->parent = parent_node; //give previous node as parent node
            }
            else { // If branch there resolve by insetting a zero length branch
                node *extra_node = new node;   // the extra node for 0 branch length
                extra_node->right = present_node->right; // the new node will get the branches of the present node
                extra_node->left = present_node->left;   // point the left daughter to be a descendent of extra node
                extra_node->right->parent = extra_node;  // point the right daughter to be a descendent of extra node
                extra_node->left->parent = extra_node;
                extra_node->parent = present_node; // it is connected to the present node
                present_node->left = extra_node; // connect to the extra node to the left
                present_node->right = new node; //create new right node
                parent_node = present_node; //store present node location
                present_node = present_node->right; //go to newly created node
                present_node->parent = parent_node; //give previous node as parent node
            }
        }
        //if right bracket move back to parent node
        else if ( temp == ')' ) {
	    read_mode = 'l';
            if (!label_string.empty() || !extra_label_string.empty()) {
		if (!label_translation.empty()) {
		    map<string,string>::const_iterator translation = label_translation.find(label_string);
		    if (translation != label_translation.end()) label_string = translation->second;
		}
    		present_node->nodelabel = nodelabels.add_string(label_string + extra_label_string);
    		label_string.clear();
		extra_label_string.clear();
		#ifdef DEBUG
		cerr << *present_node->nodelabel << endl;
		#endif //DEBUG
            }
	    if (!branch_length_annotation.empty()) {
		present_node->branchlength = atof(branch_length_annotation.c_str());
		branch_length_annotation.clear();
	    }
            present_node = present_node->parent; //move to parent node
        }
        //if colon read branch length
        else if ( temp == ':') read_mode = 'b';
        else if ( temp == ';' ) {
            if (present_node == root ) { //check so we are back at root, if so we happily finish
		if (!label_string.empty() || !extra_label_string.empty()) { // if a lable has been read
		    if (!label_translation.empty()) { // if translating labels (e.g. nexus with translations)
			map<string,string>::const_iterator translation = label_translation.find(label_string);
			if (translation != label_translation.end()) label_string = translation->second;
		    }
		    present_node->nodelabel = nodelabels.add_string(label_string + extra_label_string);
		    #ifdef DEBUG
		    cerr << *present_node->nodelabel << endl;
		    #endif //DEBUG
		}
		if (!branch_length_annotation.empty()) {
		    present_node->branchlength = atof(branch_length_annotation.c_str());
		    branch_length_annotation.clear();
		}
		if (present_node->right == 0) {
		    if (present_node->left != 0 &&  present_node->nodelabel !=0) { // Rooted on leaf
			present_node->right = new node;
			present_node->right->nodelabel=present_node->nodelabel;
			present_node->nodelabel = 0;
			present_node->right->branchlength = present_node->branchlength;
			present_node->branchlength = 0.0;
			present_node->right->parent = present_node;
		    }
		    else cerr << "Only one child to root, this may cause problems for some functions." << endl;
		}
		break;
	    }
            else { //if we are not back at the root something is wrong
                std::cerr << "Right and left brackets not balanced in tree file, removing tree." << endl; //print error massage
                destroy_tree(root); //destroy what we have built
                break; //finish this shit
            }
        }
	else if (read_mode == 'b') {
	    branch_length_annotation += temp;
        }
        //if we read semicolon we have reached the end of the tree
        else if (read_mode=='l') {
	    if (label_string.empty() && (temp == '\'' || temp == '"')) {
		char marker = temp;
		//unsigned int marker_count(1);
		infile >> temp;
		while (infile) {
		    if (temp == marker) break;
		    else label_string += temp;
		    infile >> temp;
		}
	    }
            else label_string += temp;
        }
        infile >> temp; //read next character for next loop
    }
}

int tree::print_newick_subtree( ostream& output, node *leaf, int n, bool include_br_length, bool int_node_label, map<string,string>& translate_taxa ) {
    if (leaf->left != 0) {
        output << "(";
        n=print_newick_subtree(output, leaf->left, n, include_br_length, int_node_label, translate_taxa);
    }
    else { 
        if (n<0 && leaf->nodelabel != 0) output << *leaf->nodelabel;
        else if (translate_taxa.empty()) { 
            output << n;
	    if (leaf->nodelabel!=0 ) {
		size_t comment_pos = leaf->nodelabel->find("[");
		if (comment_pos < leaf->nodelabel->length())
		    output << leaf->nodelabel->substr(leaf->nodelabel->find("["));
	    }
            ++n;
        }
	else if (leaf->nodelabel!=0 ) {
	    #ifdef DEBUG
	    cerr << "Trying to translate " << *leaf->nodelabel << endl;
	    #endif //DEBUG
	    string name;
	    size_t comment_pos = leaf->nodelabel->find("[");
	    if (comment_pos < leaf->nodelabel->length())
		name = leaf->nodelabel->substr(0,comment_pos);
	    else name = *leaf->nodelabel;
	    //cerr << name << endl;
	    for (map<string,string>::iterator i = translate_taxa.begin(); i!=translate_taxa.end(); ++i)
		if (!i->second.compare(name)) { output << i->first; break; }
	    if (comment_pos < leaf->nodelabel->length())
		output << leaf->nodelabel->substr(leaf->nodelabel->find("["));
	}
    }
    if (leaf->right != 0) {
        output << ",";
        n=print_newick_subtree(output, leaf->right, n, include_br_length, int_node_label, translate_taxa);
        output << ")";
        if (int_node_label && leaf->nodelabel !=0) output << *leaf->nodelabel;
    }
    if (leaf != root && include_br_length) output << ":" << fixed << leaf->branchlength;
    return n;
}

void tree::print_nexus_subtree( ostream& output, node *leaf, bool include_br_length, bool int_node_label ) {
    map<string,string> translate_taxa;
    print_nexus_tree_head( output, leaf, translate_taxa );
    string name = "tree1";
    print_tree_to_nexus_stream(output,leaf,name,include_br_length,int_node_label, translate_taxa);
    output << "End;" << endl;
}

void tree::print_nexus_tree_head( ostream& output, node *leaf, map<string,string>& translate_taxa ) {
    output << "#NEXUS" << endl;
    output << "Begin taxa;" << endl;
    output << "\tDimensions ntax=" << n_sub_tips( leaf ) << ";" << endl;
    output << "\tTaxlabels" << endl;
    if (translate_taxa.empty()) {
	vector<string*> tips;
	tip_names( leaf, tips);
	for (vector<string*>::iterator i=tips.begin(); i != tips.end(); ++i) {
	    size_t comment_start = (*i)->find("[");
	    string name = (*i)->substr(0,comment_start);
	    output << "\t\t" << name << endl;
	}
	output << ";" << endl << "End;" << endl;
	output << endl;
	output << "Begin trees;" << endl << "\tTranslate" << endl;
	int count = 0;
	for (vector<string*>::iterator i=tips.begin(); i != tips.end(); ++i) {
	    size_t comment_start = (*i)->find("[");
	    string name = (*i)->substr(0,comment_start);
	    output << "\t\t" << ++count << ' ' << name;
	    if (i+1 != tips.end()) output << ',' << endl;
	    stringstream ss;
	    ss << count;
	    string number = ss.str();
	    translate_taxa[number]=name;
	}
    }
    else {
	for (map<string,string>::const_iterator i=translate_taxa.begin(); i != translate_taxa.end(); ++i) {
	    output << "\t\t" << i->second << endl;
	}
	output << ";" << endl << "End;" << endl;
	output << endl;
	output << "Begin trees;" << endl << "\tTranslate" << endl;
	int count = translate_taxa.size();
	for (map<string,string>::const_iterator i=translate_taxa.begin(); i != translate_taxa.end(); ++i) {
	    output << "\t\t" << i->first << ' ' << i->second;
	    --count;
	    if (count>0)  output << ',' << endl;
	}
    }
    output << endl;
    output << ";" << endl;
}
void tree::drop_tips_on_short_branches (node* leaf, const float cut_off) {
    if (leaf->left == 0 && leaf->right == 0 && leaf->parent != 0 && leaf->branchlength < cut_off)
	drop_tips(leaf->parent,*leaf->nodelabel);
    else {
	if (leaf->left != 0) drop_tips_on_short_branches(leaf->left, cut_off);
	if (leaf->right != 0) drop_tips_on_short_branches(leaf->right, cut_off);
    }
}
char tree::drop_tips (node *leaf, const string taxa) {
    if (leaf->left == 0 && leaf->right == 0) {
        int length = taxa.length();
        string temp;
        for (int i=0; i <= length; ++i) {
            if (taxa[i] == ',' || i == length) {
                if (!temp.compare(*leaf->nodelabel)) {
                    return 'd';
                }
                temp.clear();
            }
            else temp += taxa[i];
        }
        return 'k';
    }
    else {
        char left = 'd';
        char right = 'd';
        if ( leaf->left != 0 ) { left = drop_tips( leaf->left, taxa ); }
        if ( leaf->right != 0 ) { right = drop_tips( leaf->right, taxa ); }
        if ( left == right ) return left;
        else {
            if (leaf != root) {
                if ( leaf->parent->left == leaf && right == 'k' ) {
                    leaf->parent->left = leaf->right;
                    leaf->right->parent = leaf->parent;
                    leaf->right->branchlength += leaf->branchlength;
                    leaf->right = 0;
                }
                else if ( leaf->parent->left == leaf && left == 'k' ) {
                    leaf->parent->left = leaf->left;
                    leaf->left->parent = leaf->parent;
                    leaf->left->branchlength += leaf->branchlength;
                    leaf->left = 0;
                }
                else if ( leaf->parent->right == leaf && right == 'k' ) {
                    leaf->parent->right = leaf->right;
                    leaf->right->parent = leaf->parent;
                    leaf->right->branchlength += leaf->branchlength;
                    leaf->right = 0;
                }
                else if ( leaf->parent->right == leaf && left == 'k' ) {
                    leaf->parent->right = leaf->left;
                    leaf->left->parent = leaf->parent;
                    leaf->left->branchlength += leaf->branchlength;
                    leaf->left = 0;
                }
                destroy_tree (leaf);
                return 'k';
            }
            else {
                if (left == 'k') {
                    root = leaf->left;
                    leaf->left = 0;
                }
                else if (right == 'k') {
                    root = leaf->right;
                    leaf->right = 0;
                }
                destroy_tree (leaf);
                return 'k';
            }
        }
    }
}

double tree::longest_to_tip (node *leaf) {
    double left=0;
    double right=0;
    if ( leaf == 0 ) return 0;
    if ( leaf->left != 0 ) {
        left = longest_to_tip (leaf->left);
    }
    if ( leaf->right != 0 ) {
        right = longest_to_tip (leaf->right);
    }
    if ( left > right ) { return leaf->branchlength + left; }
    else { return leaf->branchlength + right; }
}

void tree::tip_furthest_from_root (node *leaf, node_and_distance* tip_distance) {
    if ( leaf == 0 ) {
        tip_distance->tip=0;
        tip_distance->distance=0;
        return;
    }
    if ( leaf->left == 0 && leaf->right == 0) {
        tip_distance->tip = leaf;
        tip_distance->distance=leaf->branchlength;
        return;
    }
    node_and_distance left/*= new node_and_distance*/;
    node_and_distance right/*= new node_and_distance*/;
    if ( leaf->left != 0 ) {
        tip_furthest_from_root (leaf->left,&left);
    }
    if ( leaf->right != 0 ) {
        tip_furthest_from_root (leaf->right,&right);
    }
    if ( left.distance > right.distance ) {
        tip_distance->distance = leaf->branchlength + left.distance;
        tip_distance->tip = left.tip;
    }
    else {
        tip_distance->distance = leaf->branchlength + right.distance; 
        tip_distance->tip = right.tip;
    }
    //delete left;
    //delete right;
}

double tree::shortest_to_tip (node *leaf) {
    if ( leaf == 0 ) return 0;
    double left=0;
    double right=0;
    if ( leaf->left != 0 ) {
        left = shortest_to_tip (leaf->left);
    }
    if ( leaf->right != 0 ) {
        right = shortest_to_tip (leaf->right);
    }
    if ( left <= right ) { return leaf->branchlength + left; }
    else { return leaf->branchlength + right; }
}


string tree::tip_names ( node *leaf, const string& separator ) {
    if ( leaf->left == 0 && leaf->right == 0 ) {
        return *leaf->nodelabel;
    }
    else {
        string tips;
        if ( leaf->left !=0 ) tips = tip_names( leaf->left, separator );
        if ( leaf->right !=0 ) {
            if (!tips.empty()) tips += separator;
            tips += tip_names( leaf->right, separator );
        }
        return tips;
    }
}

void tree::tip_names ( node *leaf, set<string*>& tips ) {
    if ( leaf->left == 0 && leaf->right == 0 ) {
        tips.insert(leaf->nodelabel);
    }
    else {
        if ( leaf->left !=0 ) tip_names( leaf->left, tips );
        if ( leaf->right !=0 ) tip_names( leaf->right, tips );
    }
}

void tree::tip_names ( node *leaf, vector<string*>& tips ) {
    if ( leaf->left == 0 && leaf->right == 0 ) {
        tips.push_back(leaf->nodelabel);
    }
    else {
        if ( leaf->left !=0 ) tip_names( leaf->left, tips );
        if ( leaf->right !=0 ) tip_names( leaf->right, tips );
    }
}

bool tree::tip_is_in_list (const string* tip, const string* list, const char separator) {
    string test;
    int length = list->length();
    for (int i=0; i<=length; ++i) {
        if ( i  == length || list->at(i) == separator ) {
            if (!tip->compare(test)) return true;
            else test.clear();
        }
        else test+=list->at(i);
    }
    return false;
}
int tree::print_tips ( ostream& output, node *leaf, string leading_n, int n, string leading, string trailing ) {
    if ( leaf == 0 ) return n;
    if ( leaf->left != 0 || leaf->right != 0 ) {
        if ( leaf->left != 0) n=print_tips( output, leaf->left, leading_n, n, leading, trailing );
        if ( leaf->right != 0) n=print_tips( output, leaf->right, leading_n, n, leading, trailing );
    }
    else {
        if (n < 0 && leaf->nodelabel!=0) output << leading << *leaf->nodelabel << trailing;
        else {
            if (leaf->nodelabel!=0) output << leading_n << n << leading << *leaf->nodelabel << trailing;
            ++n;
        }
    }
    return n;
}

void tree::print_non_descendants( ostream& output, node *leaf, node* done, const string leading, const string trailing ) {
// the node that is given as done will be excluded from printing
    if ( leaf != 0 ) {
	if (leaf->left != 0 && leaf->left != done ) print_non_descendants( output, leaf->left, leaf, leading, trailing );
	if (leaf->right != 0 && leaf->right != done ) print_non_descendants( output, leaf->right, leaf, leading, trailing );
	if (leaf->parent != 0 && leaf->parent != done ) print_non_descendants( output, leaf->parent, leaf, leading, trailing );
	if (leaf->left == 0 && leaf->right == 0 && leaf->nodelabel != 0) {
	    output << leading << *leaf->nodelabel << trailing;
	}
    }
}
void tree::midpoint_root() {
    node* midpoint_node = root;
    double branch_length1 = find_midpoint_node(midpoint_node); // find node with midpoint
    re_root(midpoint_node); // root on that node
    // Adjust branchlengths from the root
    double root_branch = root->left->branchlength + root->right->branchlength;
    midpoint_node->branchlength = branch_length1;
    #ifdef DEBUG
    cerr << "Branch length: " << branch_length1 << " Root branch: " << root_branch <<endl << root << " vs " << midpoint_node << endl;
    #endif //DEBUG
    if (root->left == midpoint_node)
	root->right->branchlength = root_branch-branch_length1;
    else if (root->right == midpoint_node)
	root->left->branchlength = root_branch-branch_length1;
    else cerr << "Warning! Failure in midpoint rooting." << endl;

}

double tree::find_midpoint_node (node*& leaf ) {
    midpoint_data midpoint_path;
    find_longest_path (leaf,midpoint_path);
    if (midpoint_path.left == 0 || midpoint_path.right == 0 || midpoint_path.mrca == 0) return 0.0;
    #ifdef DEBUG
    cerr << "longest path is between " << *midpoint_path.left->nodelabel << " (" << midpoint_path.left_path << ") and " << *midpoint_path.right->nodelabel << " (" << midpoint_path.right_path << ")" << endl;
    #endif //DEBUG
    if (midpoint_path.left_path > midpoint_path.right_path)
	return find_mid_path(leaf,midpoint_path.left,midpoint_path.mrca,(midpoint_path.left_path+midpoint_path.right_path)/2.0);
    else
	return find_mid_path(leaf,midpoint_path.right,midpoint_path.mrca,(midpoint_path.left_path+midpoint_path.right_path)/2.0);
}

double tree::find_mid_path ( node*& midpoint_node, node* start, const node* stop, const double path_length ) {
    if (start == 0) return 0.0;
    #ifdef DEBUG
    cerr << path_length << " left to midpoint" << " (" << start << ")" << endl;
    #endif //DEBUG
    if (start->branchlength >= path_length) {
	midpoint_node = start;
	//if (start->parent == stop)
      	return path_length;
	//return start->branchlength - path_length;
    }
    else if (start->parent != 0 && start->parent != stop)
	return find_mid_path(midpoint_node, start->parent, stop, path_length-start->branchlength);
    else {
	cerr << "Warning! Could not find mid point." << endl;
	return 0.0;
    }
}

double tree::find_longest_path ( node* leaf, midpoint_data& path_data ) {
    if (leaf == 0) return 0.0;
    if (leaf->left == 0 && leaf->right == 0 && leaf->parent != 0) {
	if (leaf == leaf->parent->left) {
	    path_data.left = leaf;
	    //path_data.left_path = leaf->branchlength;
	}
	else if (leaf == leaf->parent->right) {
	    path_data.right = leaf;
	    //path_data.right_path = leaf->branchlength;
	}
	return leaf->branchlength;
    }
    else {
	midpoint_data data_left;
	double longest_left(0.0);
	midpoint_data data_right;
	double longest_right(0.0);
	if (leaf->left != 0) longest_left = find_longest_path(leaf->left,data_left);
	if (leaf->right != 0) longest_right = find_longest_path(leaf->right,data_right);
	if (longest_left+longest_right > data_left.right_path+data_left.left_path && longest_left+longest_right > data_right.right_path+data_right.left_path) { // if longest path through this node
	    path_data.mrca = leaf;
	    if (data_left.left == 0)
		path_data.left = data_left.right;
	    else if (data_left.right == 0)
		path_data.left = data_left.left;
	    else if (data_left.left_path > data_left.right_path)
		path_data.left = data_left.left;
	    else
		path_data.left = data_left.right;
	    path_data.left_path = longest_left;
	    if (data_right.left == 0)
		path_data.right = data_right.right;
            else if (data_right.right == 0)
		path_data.right = data_right.left;
	    else if (data_right.left_path > data_right.right_path)
		path_data.right = data_right.left;
	    else
		path_data.right = data_right.right;
	    path_data.right_path = longest_right;
	    //if (path_data.left_path > path_data.right_path)
	//	return path_data.left_path+leaf->branchlength;
	 //   else
	//	return path_data.right_path+leaf->branchlength;
	}
	else if (data_left.right_path+data_left.left_path > data_right.right_path+data_right.left_path) {
	    path_data = data_left;
	}
	else {
	    path_data = data_right;
	}
	// If the longest path is not through this node deside which path is the longest through this node and return it
	if (longest_left > longest_right)
	    return longest_left+leaf->branchlength;
	else
	    return longest_right+leaf->branchlength;
    }
}

/*

void tree::midpoint_root ( ) {
    node *present = find_midpoint_node( root );
    double length_left;
    double length_right;
    const int precision = 100000; // precission for seting midpoint root
    re_root(present);
    length_left = longest_to_tip (root->left);
    length_right = longest_to_tip (root->right);
    if (fabs(length_left-length_right) > (root->left->branchlength+root->right->branchlength)+(1.0/precision)) {
        std::cerr << "WARNING!!! Midpoint rooting failed!!!" << endl;
    }
    else {
        if (length_left > length_right) {
            root->left->branchlength -= (length_left-length_right)/2;
            root->right->branchlength += (length_left-length_right)/2;
        }
        else {
            root->right->branchlength -= (length_right-length_left)/2;
            root->left->branchlength += (length_right-length_left)/2;
        }
    }
}

tree::node* tree::find_midpoint_node ( node* leaf ) {
    node *present; // Pointer to move around the tree
    present = leaf; // start att node given in function call
    double goal = 0.0; // the length from tip goal
    double length_left; // the longest length to tip to the left
    double length_right; // same for right
    const int precision = 100000; // precission for seting midpoint root
    while (1) {
        length_left = longest_to_tip (present->left); // get lengths
        length_right = longest_to_tip (present->right);
        if ( (length_left+length_right)/2 > goal ) { goal = (length_left+length_right)/2; } // set the goal if higher than previous goal
        if ( present->right == 0 || present->left == 0 ) { // if at a tip
            break; // break
        }
        else if ( length_left-length_right < goal/precision && length_right-length_left < goal/precision ) { // if close to equal longest length to tip
            break;
        }
        else if ( length_left-length_right > goal/precision && length_left > goal ) { // if it is longer to tip to the left move to the left
            present = present->left;
        }
        else if ( length_right-length_left > goal/precision && length_right > goal ) { // same but right
            present = present->right;
        }
        else { // if neither the right or left path is longer we are at the midpoint
            break;
        }
    }
    return present;
}
*/
void tree::re_root ( node *leaf ) {
    if ( leaf != root && leaf != 0 && leaf->parent != 0 && leaf->parent!=root) {
        if (leaf == leaf->parent->left) {
            turn_nodes(leaf->parent,0);
            root->left = leaf;
            root->right = leaf->parent;
        }
        else {
            turn_nodes(leaf->parent,1);
            root->left = leaf->parent;
            root->right = leaf;
        }
        leaf->branchlength /= 2;
        leaf->parent->branchlength = leaf->branchlength;
        leaf->parent->parent = root;
        leaf->parent = root;
    }
}
//Funktion to assist re-rooting
void tree::turn_nodes ( node *leaf, bool clockwise ) {
    if (leaf->parent != root && leaf->parent !=0) { // if not at the root
        if (leaf == leaf->parent->left) turn_nodes (leaf->parent,0); // the left node should be turned anticlockwise
        else turn_nodes (leaf->parent,1); // the right node should be turned clockwise
        if (clockwise) turn_node_clockwise ( leaf );
        else turn_node_anticlockwise ( leaf );
    }
    else {
        if (clockwise) {
            node *temp = new node;
            temp->left = leaf->left;
            if (leaf == leaf->parent->left) { 
                leaf->left = leaf->parent->right;
                leaf->parent->right->parent=leaf;
                leaf->parent->right->branchlength += leaf->branchlength; // Since we are at the root get both branches
            }
            else {
                leaf->left = leaf->parent->left;
                leaf->parent->left->parent = leaf;
                leaf->parent->left->branchlength += leaf->branchlength;
            }
            temp->right = leaf->right;
            leaf->right = temp->left;
            leaf->parent = temp->right;
            leaf->branchlength = temp->right->branchlength;
            if (temp->right->left==0 && temp->right->right==0) leaf->nodelabel=0;
            else leaf->nodelabel = temp->right->nodelabel; //test
            delete temp;
        }
        else {
            node *temp = new node;
            temp->left=leaf->left;
            leaf->left = leaf->right;
            temp->right = leaf->right;
            if (leaf == leaf->parent->left) { 
                leaf->right = leaf->parent->right;
                leaf->parent->right->parent=leaf;
                leaf->parent->right->branchlength += leaf->branchlength;
            }
            else {
                leaf->right = leaf->parent->left;
                leaf->parent->left->parent = leaf;
                leaf->parent->left->branchlength += leaf->branchlength;
            }
            leaf->parent = temp->left;
            leaf->branchlength = temp->left->branchlength;
            if (temp->left->left==0 && temp->left->right==0) leaf->nodelabel=0;
            else leaf->nodelabel = temp->left->nodelabel; //test
            delete temp;
        }
    }
}

inline void tree::turn_node_clockwise ( node *leaf ) {
    node *temp = new node;
    temp->left = leaf->left;
    leaf->left = leaf->parent;
    temp->right = leaf->right;
    leaf->right = temp->left;
    leaf->branchlength = temp->right->branchlength;
    if (temp->right->left==0 && temp->right->right==0) leaf->nodelabel=0;
    else leaf->nodelabel = temp->right->nodelabel;
    leaf->parent = temp->right;
    delete temp;
}

inline void tree::turn_node_anticlockwise ( node *leaf ) {
    node *temp = new node;
    temp->left=leaf->left;
    leaf->left = leaf->right;
    temp->right = leaf->right;
    leaf->right = leaf->parent;
    leaf->branchlength = temp->left->branchlength;
    if (temp->left->left==0 && temp->left->right==0) leaf->nodelabel=0;
    else leaf->nodelabel = temp->left->nodelabel;
    leaf->parent = temp->left;
    delete temp;
}

string tree::not_present ( node *leaf, const string taxa ) {
    string temp;
    string not_present; 
    string node_lables = tip_names( leaf, "," );
    int length_node_lables = node_lables.length();
    int i; 
    int unsigned j;
    for (i=0; i <= length_node_lables; ++i) {
        if (node_lables[i] == ',' or i == length_node_lables) {
           char add = 'y';
           string comp_temp;
           for (j=0; j <= taxa.length(); ++j) {
                if (taxa[j] == ',' or j == taxa.length()) {
                    if (!comp_temp.compare(temp)) {
                        add = 'n';
                        break;
                    }
                    comp_temp.clear();
                }
                else comp_temp += taxa[j];
            }
            if (add =='y') {
                if (not_present.empty()) not_present = temp;
                else {
                    not_present += ',';
                    not_present += temp;
                }
            }
            temp.clear();
        }
        else temp += node_lables[i];
    }
    return not_present;
}

inline bool tree::taxon_present( node *leaf, const string taxon ) {
    string node_lables = tip_names( leaf, "," );
    int length=node_lables.length();
    string comp_taxon;
    for (int i=0; i <= length; ++i) {
        if (node_lables[i] == ',' || i==length) {
            if (!comp_taxon.compare(taxon)) return 1;
            comp_taxon.clear();
        }
        else comp_taxon += node_lables[i];
    }
    return 0;
}

void tree::change_tip_name( node *leaf, const string tip_name, const string new_tip_name) {
    if (leaf->left == 0 && leaf->right == 0 && !leaf->nodelabel->compare(tip_name)) leaf->nodelabel = nodelabels.add_string(new_tip_name);
    else {
        if (leaf->left != 0) change_tip_name ( leaf->left, tip_name, new_tip_name);
        if (leaf->right != 0) change_tip_name ( leaf->right, tip_name, new_tip_name);
    }
}

void tree::change_tip_names( const string tip_pairs ) {
    int length = tip_pairs.length();
    char mode = 'p';
    string name;
    string new_name;
    for (int i=0; i <= length; ++i) {
        if (i == length || tip_pairs[i] == ',') {
            change_tip_name( root, name, new_name );
            new_name.clear();
            name.clear();
            mode = 'p';
        }
        else if (tip_pairs[i] == '|') mode = 'n';
        else if (mode == 'n') new_name += tip_pairs[i];
        else if (mode == 'p') name += tip_pairs[i];
    }
}

void tree::multiply_br_length_subtree ( node *leaf, const float multiplier ) {
    if (leaf == 0) return;
    else {
        leaf->branchlength *= multiplier;
        if (leaf->left != 0) multiply_br_length_subtree ( leaf->left, multiplier );
        if (leaf->right != 0) multiply_br_length_subtree ( leaf->right, multiplier );
    }
}

void tree::multiply_br_length_skyline (node *leaf, rate_model& rate_mod , const float distance_to_root) {
    const double branch_end = distance_to_root+leaf->branchlength;
    leaf->branchlength = rate_mod.get_time_mod_branch_length(leaf->branchlength,distance_to_root);
    if (leaf->left != 0) multiply_br_length_skyline(leaf->left, rate_mod, branch_end);
    if (leaf->right != 0) multiply_br_length_skyline(leaf->right, rate_mod, branch_end);
}

void tree::multiply_br_length_cut_off_subtree (node *leaf, const float cut_off, const float multiplier ) {
    if (leaf == 0) return;
    else {
	float next_cut_off = cut_off - leaf->branchlength;
	if (leaf->branchlength < cut_off) leaf->branchlength *= multiplier;
	else if ( cut_off > 0.0 ) leaf->branchlength = cut_off * multiplier + leaf->branchlength - cut_off;
	if (leaf->left != 0) multiply_br_length_cut_off_subtree (leaf->left, next_cut_off, multiplier);
	if (leaf->right != 0) multiply_br_length_cut_off_subtree (leaf->right, next_cut_off, multiplier);
    }
}
void tree::multiply_br_length_clades ( const vector<string> &taxon_sets, rate_model& rate_mod) {
    unsigned int order(0);
    for (vector<string>::const_iterator i = taxon_sets.begin(); i != taxon_sets.end(); ++i) {
	clades.insert(pair<node*,unsigned int>(most_recent_common_ancestor(*i),order));
	++order;
    }
    multiply_br_length_clades( root, rate_mod,  rate_mod.get_n_rates());
};

void tree::multiply_br_length_clades ( node* leaf, rate_model& rate_mod, unsigned int present_clade) {
    leaf->branchlength = leaf->branchlength*rate_mod.get_rate_clade(present_clade);
    map<node*,unsigned int>::const_iterator present_node = clades.find(leaf);
    if (present_node != clades.end()) 
	present_clade = present_node->second;
    if (leaf->left != 0) multiply_br_length_clades( leaf->left, rate_mod, present_clade );
    if (leaf->right != 0) multiply_br_length_clades( leaf->right, rate_mod, present_clade );
}

void tree::set_br_length_subtree ( node *leaf, const float value ) {
    if (leaf == 0) return;
    else {
        leaf->branchlength = value;
        if (leaf->left != 0) set_br_length_subtree ( leaf->left, value );
        if (leaf->right != 0) set_br_length_subtree ( leaf->right, value );
    }
}

double tree::sum_br_length( node *leaf ) {
    if (leaf == 0 ) return 0;
    return sum_br_length(leaf->left) + sum_br_length(leaf->right) + leaf->branchlength;
    /*else if (leaf->left != 0 && leaf->right != 0) return sum_br_length(leaf->left) + sum_br_length(leaf->right) + leaf->branchlength;
    else if (leaf->left != 0 && leaf->right == 0) return sum_br_length(leaf->left) + leaf->branchlength;
    else if (leaf->left == 0 && leaf->right != 0) return sum_br_length(leaf->right) + leaf->branchlength;
    */
}

void tree::print_distances_to_all_for_each_taxa ( const string value_divider, const string taxon_divider ) {
    string taxon_string = tip_names(",");
    int length = taxon_string.length();
    string taxon1;
    for (int i=0; i <= length; ++i) {
        if (taxon_string[i] == ',' || i == length) {
            string taxon2;
            double distance=0;
            for (int j=0; j <= length; ++j) {
                if (taxon_string[j] == ',' || j == length) {
                    if (taxon1.compare(taxon2)) {
                        distance += two_taxa_distance( nodelabels.find_string(taxon1), nodelabels.find_string(taxon2) );
                    }
                    taxon2.clear();
                }
                else taxon2 += taxon_string[j];
            }
            std::cout << taxon1 << value_divider << distance << taxon_divider;
            taxon1.clear();
        }
        else taxon1 += taxon_string[i];
    }
}

double tree::two_taxa_distance ( const string* taxon1, const string* taxon2 ) {
    node* tip1 = find_taxon_tip (root, taxon1);
    node* tip2 = find_taxon_tip (root, taxon2);
    node* ancestor1 = tip1->parent; 
    node* ancestor2 = tip2->parent;
    double distance1 = tip1->branchlength; 
    double distance2 = tip2->branchlength;
    do {
        do {
            if (ancestor1 != ancestor2 && ancestor2 != root && ancestor2 != 0) {
                distance2 += ancestor2->branchlength;
                ancestor2 = ancestor2->parent;
            }
            if (ancestor1 == ancestor2) return distance1 + distance2;
        } while (ancestor2 != root && ancestor2 != 0);
        ancestor2 = tip2->parent;
        distance2 = tip2->branchlength;
        if (ancestor1 != root && ancestor1 != 0) {
            distance1 += ancestor1->branchlength;
            ancestor1 = ancestor1->parent;
        }
    } while ( ancestor1 != 0 && ancestor2 != 0); // A bit risky but they should both go to the root at least
    return 0;
}

tree::node * tree::most_recent_common_ancestor ( const string* taxon1, const string* taxon2 ) {
    vector<node*> nodes;
    nodes.push_back(find_taxon_tip (root, taxon1));
    nodes.push_back(find_taxon_tip (root, taxon2));
    node* mrca =  most_recent_common_ancestor ( nodes );
    return mrca;
}

tree::node * tree::most_recent_common_ancestor ( const string& taxa ) {
    string taxon;
    int str_length=taxa.length();
    vector<node*> nodes;
    for (int i=0; i<=str_length; ++i) {
        if (taxa[i]==',' || i == str_length) {
	    nodes.push_back(find_taxon_tip (root, taxon));
            taxon.clear();
        }
        else taxon += taxa[i];
    }
    node* mrca =  most_recent_common_ancestor ( nodes );
    return mrca;
}

tree::node * tree::most_recent_common_ancestor ( set<string*> taxa ) {
    vector<node*> nodes;
    #ifdef DEBUG
    cerr << "N taxa: " << taxa.size() << endl;
    #endif //DEBUG
    for (set<string*>::iterator i=taxa.begin(); i != taxa.end(); ++i) {
	nodes.push_back(find_taxon_tip (root, *i));
    }
    return most_recent_common_ancestor ( nodes );
}

tree::node * tree::most_recent_common_ancestor ( vector<node*>& nodes ) {
    if (nodes.empty())  return 0;
    else if (nodes.size()==1) return nodes[0];
    node* mrca = nodes[0];
    node* tip1;
    node* tip2;
    for (vector<node*>::const_iterator i=nodes.begin()+1; i!=nodes.end(); ++i) { // go through nodes
	tip1=mrca; // tip1 one is the deepest node found so far
	if (*i != 0) tip2=*i; // tip2 is the node in order
	else continue;
        do { // for each step back in the ancestral lineage of one, check if that ancestor is in lineage of tip 2
            do { 
                if (tip1 != tip2 && tip2 != root && tip2 != 0) tip2 = tip2->parent;
                if (tip1 == tip2) { mrca=tip1; break; }
            } while (tip2 != root && tip2 != 0);
            if ( tip1 == tip2 ) break;
            if (tip1 != root && tip1 != 0) {
                tip1 = tip1->parent;
                tip2=*i;
            }
        } while (tip1 != tip2 && tip1 != 0);
    }
    return mrca;
}

bool tree::add_clades_not_in_set (node* leaf, set<node*>& nodes, set<node*>& clades) {
/* add the deepest node that is not an ancestor to any node in nodes, to clades. If leaf is deepest node it is not added, so if no nodes are added to clades either non of the nodes in nodes were present (return true) or all descendants of leaf were part of nodes*/
    if (leaf == 0) return true;
    bool left(true);
    bool right(true);
    if (leaf->left != 0) left = add_clades_not_in_set (leaf->left, nodes, clades);
    if (leaf->right != 0) right = add_clades_not_in_set (leaf->right, nodes, clades);
    if (left && !right) { clades.insert(leaf->left); }
    else if (right && !left) { clades.insert(leaf->right); }
    else {
	if (nodes.find(leaf) != nodes.end()) { return false; }
    }
    return left && right;
}

tree::node * tree::clade_stem_based ( vector<node*>& ingroup, set<node*>& outgroup ) {
/* Returns the end node of the branch that exclude all in outgroup but include all in ingroup, returns 0 if incompatable */
    if (ingroup.empty() || outgroup.empty())  return 0;
    set<node*> nested;
    if (add_clades_not_in_set(root,outgroup,nested)) {
	if (nested.empty()) return root;
    }
    node* mrca = 0;
    for (vector<node*>::const_iterator i=ingroup.begin(); i!=ingroup.end(); ++i) {
	node* tip = *i;
	while (tip != root && tip != 0) {
	    if (nested.find(tip) != nested.end()) {
		if (mrca != 0 && tip != mrca) return 0;
	    }
	    mrca = tip;
	}
    }
    return mrca;
} 

string tree::tips_in_clade_stem_based( const string& taxa, const string& separator) {
    vector<node*> ingroup;
    set<node*> outgroup;
    char mode = 's';
    int str_length=taxa.length();
    vector<node*> nodes;
    string taxon;
    for (int i=0; i<=str_length; ++i) {
        if (i == str_length || taxa[i]==',' || taxa[i]==':') {
            if (mode == 's') ingroup.push_back(find_taxon_tip (root, taxon));
	    else if (mode == 'o') outgroup.insert(find_taxon_tip (root, taxon));
            taxon.clear();
	    if (taxa[i]==':') mode = 'o';
        }
        else taxon += taxa[i];
    }
    return tip_names(clade_stem_based(ingroup, outgroup),separator);
}

tree::node * tree::find_taxon_tip ( node *leaf, string taxon ) {
    if (leaf == 0) return 0;
    else if ( leaf->left == 0 && leaf->right == 0 ) {
        if (!taxon.compare(*leaf->nodelabel)) return leaf;
        else return 0;
    }
    else {
        node *return_node=0;
        if (leaf->left != 0) return_node = find_taxon_tip( leaf->left,taxon );
        if (return_node == 0 && leaf->right != 0) return_node = find_taxon_tip( leaf->right,taxon );
        return return_node;
    }
}

tree::node * tree::find_taxon_tip ( node *leaf, const string* taxon ) {
    if (leaf == 0) return 0;
    else if ( leaf->left == 0 && leaf->right == 0 ) {
        if (leaf->nodelabel == taxon) return leaf;
        else return 0;
    }
    else {
        node *return_node=0;
        if (leaf->left != 0) return_node = find_taxon_tip( leaf->left,taxon );
        if (return_node == 0 && leaf->right != 0) return_node = find_taxon_tip( leaf->right,taxon );
        return return_node;
    }
}

void tree::print_distance_to_root ( node* leaf, double distance, int n_nodes, string value_sep, string row_sep) {
    if (leaf != 0) {
        if ( leaf->left == 0 && leaf->right == 0 && leaf->nodelabel != 0 ) std::cout << *leaf->nodelabel << value_sep << n_nodes << value_sep << distance + leaf->branchlength << row_sep;
        else {
            if ( leaf->left != 0 ) print_distance_to_root ( leaf->left, distance + leaf->branchlength, n_nodes+1, value_sep, row_sep );
            if ( leaf->right != 0 ) print_distance_to_root ( leaf->right, distance + leaf->branchlength, n_nodes+1, value_sep, row_sep );
        }
    }
}

void tree::print_branch_lengths ( string separator, bool print_root ) {
    if (root->left != 0) print_branch_lengths( root->left, separator );
    if (root->right != 0) print_branch_lengths( root->right, separator );
    if (print_root) {
        std::cout << root->branchlength;
        if (!separator.compare("\n")) std::cout << endl;
    }
    if (separator.compare("\n")) std::cout << endl;
}

void tree::print_branch_lengths ( node* leaf, string separator ) {
    if (leaf != 0 ) {
        if (leaf->left != 0) print_branch_lengths( leaf->left, separator );
        if (leaf->right != 0) print_branch_lengths( leaf->right, separator );
        std::cout << leaf->branchlength << separator;
    }
}

tree::node* tree::longest_branch ( node* leaf ) {
    if (leaf->left==0 && leaf->right==0) return leaf;
    node* longest_left = 0;
    node* longest_right = 0;
    if (leaf->left != 0) longest_left = longest_branch( leaf->left);
    if (leaf->right != 0) longest_right = longest_branch( leaf->right);
    if (longest_left != 0 && longest_right != 0) {
        if (longest_left->branchlength > longest_right->branchlength && longest_left->branchlength > leaf->branchlength) return longest_left;
        else if (longest_right->branchlength > longest_left->branchlength && longest_right->branchlength > leaf->branchlength) return longest_right;
        else return leaf;
    }
    else if (longest_left != 0 ) {
        if (longest_left->branchlength > leaf->branchlength) return longest_left;
        else return leaf;
    }
    else if (longest_right != 0 ) {
        if (longest_right->branchlength > leaf->branchlength) return longest_right;
        else return leaf;
    }
    else return leaf;
}

tree::node* tree::shortest_branch ( node* leaf ) {
    if (leaf->left==0 && leaf->right==0) return leaf;
    node* shortest_left = 0;
    node* shortest_right = 0;
    if (leaf->left != 0) shortest_left = shortest_branch( leaf->left);
    if (leaf->right != 0) shortest_right = shortest_branch( leaf->right);
    if (shortest_left != 0 && shortest_right != 0) {
        if (shortest_left->branchlength < shortest_right->branchlength && shortest_left->branchlength < leaf->branchlength) return shortest_left;
        else if (shortest_right->branchlength < shortest_left->branchlength && shortest_right->branchlength < leaf->branchlength) return shortest_right;
        else return leaf;
    }
    else if (shortest_left != 0 ) {
        if (shortest_left->branchlength < leaf->branchlength) return shortest_left;
        else return leaf;
    }
    else if (shortest_right != 0 ) {
        if (shortest_right->branchlength < leaf->branchlength) return shortest_right;
        else return leaf;
    }
    else return leaf;
}

void/*tree::tree*/ tree::prune_clade ( node* leaf,tree* pruned_clade ) {
    //tree pruned_clade;
    if (leaf == root) {
        std::cerr << "Cannot prune tree at root returning tree with no tips;";
        return /*pruned_clade*/;
    }
    else {
        (*pruned_clade).~tree();
        (*pruned_clade).root = leaf;
        (*pruned_clade).root->branchlength=0;
        node* drop_node = leaf->parent;
        if (drop_node == root) {
            node* temp = root;
            if (leaf == root->left) root = root->right;
            else root = root->left;
            root->branchlength=0;
            delete temp;
        }
        else {
            node* daughter=0;
            if (drop_node->left == leaf) daughter=drop_node->right;
            else if (drop_node->right == leaf) daughter=drop_node->left;
            if (drop_node->parent->left == drop_node) drop_node->parent->left = daughter;
            else if (drop_node->parent->right == drop_node) drop_node->parent->right = daughter;
            daughter->parent = drop_node->parent;
            daughter->branchlength += drop_node->branchlength;
            delete drop_node;
        }
        //return pruned_clade;
    }
}

void/*tree::tree*/ tree::split_at_longest_branch ( node* leaf, tree* tree ) {
    node* split_node = longest_branch(leaf);
    if (split_node == root) {
        //tree return_tree;
        std::cerr << "Root has longest branch, cutting of root branch." << endl;
        root->branchlength = 0;
        return /*return_tree*/;
    }
    /*return*/ prune_clade( split_node, tree );
}
void/*tree::tree*/ tree::split_at_longest_branch_unrooted ( node* leaf, tree* tree ) {
    if (leaf->left == 0 || leaf->right == 0) {
        std::cerr << "Basal node must be bifurcating. No split made." << endl;
        return;
    }
    node* split_node = longest_branch(leaf);
    if (split_node->branchlength >= (leaf->left->branchlength + leaf->right->branchlength)) {
        if (split_node == root) { 
            //tree return_tree;
            std::cerr << "Root has longest branch, cutting of root branch. Looking for longest branch again." << endl;
            root->branchlength = 0;
            split_node = longest_branch(leaf);
            if (split_node == root) {
                std::cerr << "Root has longest branch again. Will split at root." << endl;
                int n_left=n_sub_tips( leaf->left );
                int n_right=n_sub_tips( leaf->right );
                if (n_left > n_right) split_node = leaf->right;
                else split_node = leaf->left;
            }
        }
    }
    else {
        int n_left=n_sub_tips( leaf->left );
        int n_right=n_sub_tips( leaf->right );
        if (n_left > n_right) split_node = leaf->right;
        else split_node = leaf->left;
    }
    prune_clade( split_node, tree );
}

void tree::split_at_midpoint ( node* leaf, tree* tree ) {
    if (leaf->left == 0 || leaf->right == 0) {
        std::cerr << "Basal node must be bifurcating. No split made." << endl;
        return;
    }
    node* split_node = leaf;
    find_midpoint_node ( split_node );
    if (split_node == root) {
        int n_left=n_sub_tips( leaf->left );
        int n_right=n_sub_tips( leaf->right );
        if (n_left > n_right) prune_clade( leaf->right, tree );
        else prune_clade( leaf->left, tree );
    }
    else prune_clade( split_node, tree );
}

void tree::print_svg_autoscale (bool scalebar, bool node_lable, const unsigned int font_size, const string& tip_color ) {
    std::cout << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
    print_svg_autoscale_no_head(scalebar, node_lable, font_size, tip_color, false);
    std::cout << "</svg>" << endl;
}

void tree::print_svg_autoscale_no_head (bool scalebar, bool node_lable, const unsigned int font_size, const string& tip_color, bool for_html ) {
    node_and_distance x_lim_node;
    tip_furthest_from_root (root, &x_lim_node);
    float height = n_tips()*font_size;
    if (scalebar) height += 2*font_size;
    if (height < 100) height = 100;
    float width = height/1.618;
    unsigned int stroke_width = (font_size/20)+1;
    //if (for_html) std::cout << "<svg width=\"" << width+10 <<"\" height=\"" << height << "\">" << endl;
    print_svg_no_head (scalebar, node_lable, width, height, 5.0, stroke_width, font_size, "Arial", tip_color, for_html);
    //if (for_html) std::cout << "</svg>" << endl;
}

void tree::print_svg ( bool scalebar, bool node_lable, const float width, const float height, const float offset, const unsigned int stroke_width, const int font_size, string font, const string& tip_color ) {
    std::cout << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
    print_svg_no_head( scalebar, node_lable, width, height, offset, stroke_width, font_size, font, tip_color, false);
    std::cout << "</svg>" << endl;
}

void tree::print_svg_no_head ( bool scalebar, bool node_lable, const float width, const float height, const float offset, const unsigned int stroke_width, const int font_size, string font, const string& tip_color, bool for_html ) {
    node_and_distance x_lim_node;
    tip_furthest_from_root (root, &x_lim_node);
    float x_unit = (width-5-(x_lim_node.tip->nodelabel->length()*font_size/1.6))/x_lim_node.distance;
    float y_unit;
    if (scalebar) y_unit = (height-3*font_size)/(n_tips()+1);
    else y_unit = (height-font_size)/(n_tips()+1);
    poly_line_start stem;
    #ifdef DEBUG
    cerr << "Tip colors: " << tip_color << endl;
    #endif // DEBUG
    if (for_html) std::cout << "<svg width=\"" << width+10 <<"\" height=\"" << height << "\">" << endl;
    print_svg_subtree (root, &stem, node_lable, 1, 5.0, x_unit, y_unit, offset, stroke_width, font_size, font, tip_color);
    if (scalebar) {
        int i=0;
        while ((x_lim_node.distance/5) * pow(10.0,i) <= 9.5) ++i;
        if (i == 0) while ((x_lim_node.distance/5) * pow(10.0,i) >= 99.5) --i;
        int  int_scale = x_lim_node.distance/5 * pow(10.0,i);
        if ((x_lim_node.distance/5 * pow(10.0,i))-int_scale > 0.5) ++int_scale;
        float length_scalebar = int_scale/pow(10.0,i);
        std::cout << "<line x1=\"" << (width/2)-(length_scalebar*x_unit/2)+5 << "\" y1=\"" << height-2*font_size << "\" x2=\"" << (width/2)+(length_scalebar*x_unit/2)+5
	    << "\" y2=\"" << height-2*font_size << "\" style=\"fill:none;stroke:black;stroke-width:" << stroke_width << "\" />" << endl;
        std::cout << "<text x=\"" << (width/2)-(length_scalebar*x_unit/2)+5 << "\" y=\"" << height-0.5*font_size << "\" font-family=\"" << font << "\" font-size=\""
	    << font_size << "\">" << length_scalebar << "</text>" << endl;
    }
    if (for_html) std::cout << "</svg>" << endl;
}

void tree::print_svg_subtree ( node* leaf, poly_line_start* branch_start, bool node_lable, int n, float x, const float x_unit, const float y_unit,
                               const float offset, const int stroke_width, const int font_size, string font, const string& tip_color ) {
    x += leaf->branchlength * x_unit;
    float y;
    if (leaf->left == 0 && leaf->right == 0) {
        y = n * y_unit;
        std::cout << "<text x=\"" << (x+offset)  << "\" y=\"" << y+(font_size/3.0) << "\" font-family=\"" << font << "\" font-size=\"" << font_size ;
	if (!tip_color.empty()) {
	    unsigned int length = tip_color.length();
	    string color;
	    string temp;
	    for (unsigned int i=0; i<length+1; ++i) {
                if (i == length || tip_color.at(i) == ',' || tip_color.at(i) == ')') {
                    if (!leaf->nodelabel->compare(temp)) {
                        std::cout << "\" style=\"fill:" << color << ';';
                    }
		    temp.clear();
                }
		else if (tip_color.at(i) == '(') {
		    color.clear();
		    color = temp;
		    temp.clear();
		}
		//else if (tip_color.at(i) == ')') 
		else temp += tip_color.at(i);
	    }
	}
	if (leaf->nodelabel != 0) std::cout << "\">" << *leaf->nodelabel << "</text>" << endl;
        ++n;
    }
    poly_line_start* right_branch = new poly_line_start;
    poly_line_start* left_branch = new poly_line_start;
    if (leaf->left != 0) {
        print_svg_subtree( leaf->left, left_branch, node_lable, n, x, x_unit, y_unit, offset, stroke_width, font_size, font, tip_color );
        n = left_branch->n;
    }
    if (leaf->right != 0) {
        print_svg_subtree( leaf->right, right_branch, node_lable, n, x, x_unit, y_unit, offset, stroke_width, font_size, font, tip_color );
        n = right_branch->n;
    }
    if (left_branch->line_start.compare(0,9,"<polyline") == 0 || right_branch->line_start.compare(0,9,"<polyline") == 0) { //
        if (left_branch->line_start.compare(0,9,"<polyline") == 0 and right_branch->line_start.compare(0,9,"<polyline") == 0) {
            y = left_branch->y + ((right_branch->y - left_branch->y)/2);
            std::cout << left_branch->line_start << " " << x << "," << left_branch->y << " " << x << "," << y << "\"\nstyle=\"fill:none;stroke:black;stroke-width:" << stroke_width << "\" />" << endl;
            std::cout << right_branch->line_start << " " << x << "," << right_branch->y << " " << x << "," << y << "\"\nstyle=\"fill:none;stroke:black;stroke-width:" << stroke_width << "\" />" << endl;
        }
        else if (left_branch->line_start.compare(0,9,"<polyline") == 0) {
            y = left_branch->y;
            std::cout << left_branch->line_start << " " << x << "," << y << "\"\nstyle=\"fill:none;stroke:black;stroke-width:" << stroke_width << "\" />" << endl;
        }
        else {
            y = right_branch->y;
            std::cout << right_branch->line_start << " " << x << "," << y << "\"\nstyle=\"fill:none;stroke:black;stroke-width:" << stroke_width << "\" />" << endl;
        }
    }
    if (leaf->nodelabel != 0 && (leaf->left != 0 || leaf->right != 0)) {
	cout << "<text x=\"" << (x - (leaf->branchlength/2)*x_unit - leaf->nodelabel->length()*(font_size/3.2)) << "\" y=\"" << y - stroke_width << "\" font-family=\"" << font << "\" font-size=\""
	<< font_size << "\">" << *leaf->nodelabel << "</text>" << endl;
    }
    delete right_branch;
    delete left_branch;
    branch_start->n = n;
    branch_start->y = y;
    branch_start->line_start = "<polyline points=\"";
    std::stringstream convert_x (stringstream::in | stringstream::out);
    convert_x << x;
    branch_start->line_start += convert_x.str();
    branch_start->line_start += ",";
    std::stringstream convert_y (stringstream::in | stringstream::out);
    convert_y << y;
    branch_start->line_start += convert_y.str();
}

void tree::rand_topology ( const int n  ) {
    destroy_tree(root->right);
    root->right = 0;
    destroy_tree(root->left);
    root->left = 0;
    vector<node*> tips;
    if (n<1) return;
    for (int i=1; i<=n; ++i) {
	tips.push_back(new node);
        stringstream number;
        number << i;
	tips.back()->nodelabel = nodelabels.add_string(number.str());
    }
    for (int i=n; i>2; --i) {
        int no1 = rand() % i;
        int no2 = rand() % (i-1);
        if (no2 >= no1) ++no2;
        node* new_node= new node;
	new_node->nodelabel=0;
	new_node->branchlength=0;
	tips[no1]->parent = new_node;
        new_node->left = tips[no1];
	tips[no2]->parent = new_node;
        new_node->right = tips[no2];
	tips[no2] = new_node;
	tips.erase(tips.begin()+no1);
    }
    if (tips[0] != 0 && tips.size()==2) {
	tips[0]->parent=root;
        root->left=tips[0];
        tips[1]->parent=root;
        root->right=tips[1];
    }
    else if (tips[0] != 0 && tips.size()==1) {
        root->nodelabel = tips[0]->nodelabel;
        root->left=0;
        root->right=0;
    }
    else {
	root->left=0;
        root->right=0;
	for(vector<node*>::iterator i=tips.begin();i!=tips.end();++i){
	    destroy_tree(*i);
	}
    }
}
/*
void tree::add_Yule_node_depths (node* leaf) {
    unsigned int n_tips = n_sub_tips(leaf);
    double node_depth(0.0);
    set<node*> lineages;
    if (leaf->left != 0 && leaf->left->left != 0 && leaf->left->right != 0) lineages.insert(leaf->left);
    if (leaf->right != 0 && leaf->right->left != 0 && leaf->right->right != 0) lineages.insert(leaf->right);
    for (unsigned int node_no = 2; node_no < n_tips; ++node_no) { // for all nodes
	node_depth += log(((double) rand()/(RAND_MAX))+1)/(node_no); // get random node deapth
	unsigned int rand_node(0);
       	if (lineages.size()>0) rand_node = rand() % (lineages.size()); // select random node
	double dist_from_root(0.0);
	#ifdef DEBUG
	cerr << "node_depth: " << node_depth << " Node selected: " << rand_node << endl;
	#endif
	// get randomly selected node
	node* selected_node;
	set<node*>::const_iterator i = lineages.begin();
	for (; i !=lineages.end() && rand_node > 0; ++i) --rand_node;
	selected_node = *i;
	///////
	node* parent = selected_node->parent;
	#ifdef DEBUG
	cerr << "Made it here (" << rand_node << ')' << endl;
	print_tips(cerr,selected_node);
	cerr << endl;
	#endif //DEBUG
	while (parent != 0 && parent != leaf) { // get distance from root
	    dist_from_root += parent->branchlength;
	    parent = parent->parent;
	}
	selected_node->branchlength = node_depth-dist_from_root; // assign branch length
	lineages.erase(selected_node); // remove node that has gotten branch length
	// add its daughters as possible nodes to assign branch lengths to
	if (selected_node->left != 0 && selected_node->left->left != 0 && selected_node->left->right != 0)
	    lineages.insert(selected_node->left);
	if (selected_node->right != 0 && selected_node->right->left != 0 && selected_node->right->right != 0)
	    lineages.insert(selected_node->right);
	#ifdef DEBUG
	cerr << "Number of nodes: " << lineages.size() << endl;
	#endif //DEBUG
    }
    #ifdef DEBUG
    cerr << "Will set tip labels." << endl;
    #endif //DEBUG
    set_tip_branch_lengths(leaf,node_depth+log(((double) rand()/(RAND_MAX))+1)/n_tips);
}
*/
void tree::add_rand_node_depths (node* leaf, char model) {
    unsigned int n_tips = n_sub_tips(leaf);
    double node_depth(0.0);
    set<node*> lineages;
    if (leaf->left != 0 && leaf->left->left != 0 && leaf->left->right != 0) lineages.insert(leaf->left);
    if (leaf->right != 0 && leaf->right->left != 0 && leaf->right->right != 0) lineages.insert(leaf->right);
    forward_list<double> node_depths;
    if (model == 'c') {
    node_depths.push_front(-log((double) rand()/(RAND_MAX))/(n_tips*(n_tips-1)/2.0));
      	#ifdef DEBUG
	    cerr << "node_depth: " << node_depths.front() << endl;
        #endif
	for (unsigned int i=n_tips-1; i > 1; --i) {
	    node_depths.push_front(-log((double) rand()/(RAND_MAX))/(i*(i-1)/2.0)+node_depths.front());
	    #ifdef DEBUG
		cerr << "node_depth: " << node_depths.front() << endl;
	    #endif
	}
	if (!node_depths.empty()) {
	    double root_age = node_depths.front();
	    forward_list<double>::iterator previous = node_depths.begin();
	    forward_list<double>::iterator i=node_depths.begin();
	    ++i;
	    for (; i != node_depths.end(); ++i) {
		*previous = root_age-*i;
		previous = i;
	    }
	    *previous = root_age;
	}
    }
    for (unsigned int node_no = 2; node_no < n_tips; ++node_no) { // for all nodes
	if (model == 'c') {
	    node_depth = node_depths.front();
	    node_depths.pop_front();
	}
	else {
	    node_depth += log(((double) rand()/(RAND_MAX))+1)/(node_no); // get random node deapth
	}
	unsigned int rand_node(0);
       	if (lineages.size()>0) rand_node = rand() % (lineages.size()); // select random node
	double dist_from_root(0.0);
	#ifdef DEBUG
	cerr << "node_depth: " << node_depth << " Node selected: " << rand_node << endl;
	#endif
	// get randomly selected node
	node* selected_node;
	set<node*>::const_iterator i = lineages.begin();
	for (; i !=lineages.end() && rand_node > 0; ++i) --rand_node;
	selected_node = *i;
	///////
	node* parent = selected_node->parent;
	#ifdef DEBUG
	cerr << "Made it here (" << rand_node << ')' << endl;
	print_tips(cerr,selected_node);
	cerr << endl;
	#endif //DEBUG
	while (parent != 0 && parent != leaf) { // get distance from root
	    dist_from_root += parent->branchlength;
	    parent = parent->parent;
	}
	selected_node->branchlength = node_depth-dist_from_root; // assign branch length
	lineages.erase(selected_node); // remove node that has gotten branch length
	// add its daughters as possible nodes to assign branch lengths to
	if (selected_node->left != 0 && selected_node->left->left != 0 && selected_node->left->right != 0)
	    lineages.insert(selected_node->left);
	if (selected_node->right != 0 && selected_node->right->left != 0 && selected_node->right->right != 0)
	    lineages.insert(selected_node->right);
	#ifdef DEBUG
	cerr << "Number of nodes: " << lineages.size() << endl;
	#endif //DEBUG
    }
    #ifdef DEBUG
    cerr << "Will set tip labels." << endl;
    #endif //DEBUG
    if (model == 'c') set_tip_branch_lengths(leaf,node_depths.front());
    else set_tip_branch_lengths(leaf,node_depth+log(((double) rand()/(RAND_MAX))+1)/n_tips);
}

void tree::set_tip_branch_lengths (node* leaf, const double tip_branch) {
    if (leaf == 0) return;
    if (leaf->left == 0 && leaf->right == 0 && leaf->parent != 0 ) leaf->branchlength = tip_branch;
    if (leaf->left != 0) set_tip_branch_lengths (leaf->left,tip_branch-leaf->branchlength);
    if (leaf->right != 0) set_tip_branch_lengths (leaf->right,tip_branch-leaf->branchlength);
}

tree::node* tree::sim_tree_bd ( const double b, const double d, const double max, node* parent ) {
    if (b <= 0.0 || max <= 0.0) return 0;
    node* new_node = new node;
    new_node->parent = parent;
    double t_speciation(marth::expr(b));
    double t_extinction;
    if (d > 0.0) t_extinction = marth::expr(d);
    else t_extinction = t_speciation;
    if (t_extinction < t_speciation) {
	if (t_extinction > max) new_node->branchlength = max;
	else new_node->branchlength = t_extinction;
    }
    else {
	if (t_speciation > max) new_node->branchlength = max;
	else {
	    new_node->branchlength = t_speciation;
	    new_node->left = sim_tree_bd(b,d,max-t_speciation, new_node);
	    new_node->right = sim_tree_bd(b,d,max-t_speciation, new_node);
	}
    }
    #ifdef DEBUG
    cerr << new_node->branchlength << " - " << max << endl;
    #endif //DEBUG
    return new_node;
}

unsigned int tree::set_tip_lables( node* leaf, unsigned int n ) {
    if (leaf == 0) return n;
    if (leaf->left != 0 || leaf->right != 0) {
	if (leaf->left != 0) n = set_tip_lables(leaf->left, n);
	if (leaf->right != 0) n = set_tip_lables(leaf->right, n);
    }
    else if (leaf->nodelabel==0) {
	leaf->nodelabel = nodelabels.add_string("Tip" + to_string(n));
	++n;
    }
    return n;
}

void tree::add_branches_to_tree( node* leaf, const double b, const double d ) {
    if (leaf != 0) {
	double t_speciation(marth::expr(b));
	if (t_speciation < leaf->branchlength) {
	    node* new_node = new node;
	    new_node->branchlength = t_speciation;
	    new_node->parent = leaf->parent;
	    if (new_node->parent->left == leaf) new_node->parent->left = new_node;
	    else if (new_node->parent->right == leaf) new_node->parent->right = new_node;
	    leaf->branchlength = leaf->branchlength-t_speciation;
	    leaf->parent = new_node;
	    new_node->left = leaf;
	    new_node->right = sim_tree_bd(b,d, longest_to_tip(leaf)+leaf->branchlength, new_node);
	    if (leaf == root) root = new_node;
	}
	if (leaf->left != 0) add_branches_to_tree( leaf->left, b, d );
	if (leaf->right != 0) add_branches_to_tree( leaf->right, b, d );
    }
}


/*** Sub to count the number of nodes with support value above or equal to given value ***/
int tree::n_supported (node* leaf, float cutoff) {
    if (leaf->left != 0 && leaf->right != 0) {
        if (leaf->nodelabel != 0 && atof(leaf->nodelabel->c_str())>cutoff) return 1+n_supported(leaf->left, cutoff)+n_supported(leaf->right, cutoff);
        else return n_supported(leaf->left, cutoff)+n_supported(leaf->right, cutoff);
    }
    else if (leaf->left != 0) {
        if (leaf->nodelabel != 0 && atof(leaf->nodelabel->c_str())>cutoff) return 1+n_supported(leaf->left, cutoff);
        else return n_supported(leaf->left, cutoff);
    }
    else if (leaf->right != 0) {        
        if (leaf->nodelabel != 0 && atof(leaf->nodelabel->c_str())>cutoff) return 1+n_supported(leaf->right, cutoff);
        else return n_supported(leaf->right, cutoff);
    }
    else return 0;
}

bool tree::is_monophyletic ( node* leaf, string& taxa ) {
    if ( taxa.empty() ) return false;
    set<string*> query_taxa;
    int unsigned length = taxa.length();
    string taxon;
    for (int unsigned i=0; i <= length; ++i) {
        if (taxa[i]==',' || i==length) {
	    string* query = nodelabels.find_string(taxon);
	    if (query != 0) query_taxa.insert(query);
	    taxon.clear();
	}
	else taxon += taxa[i];
    }
    return is_monophyletic( leaf, query_taxa);
}

bool tree::is_monophyletic ( node* leaf, set<string*>& taxa ) {
    if (leaf->left == 0 && leaf->right == 0) {
	if (taxa.find(leaf->nodelabel)!=taxa.end()) return true;
	else return false; 
    }
    bool mono_left;
    bool mono_right;
    if (leaf->left != 0) mono_left = is_monophyletic( leaf->left, taxa );
    if (leaf->right != 0) mono_right = is_monophyletic( leaf->right, taxa );
    if (mono_left && mono_right) return true;
    else if (mono_left && n_sub_tips(leaf->left) >= taxa.size()) return true;
    else if (mono_right && n_sub_tips(leaf->right) >= taxa.size()) return true;
    else return false;
}

char tree::is_monophyletic_unrooted (node* leaf, set<string*>& taxa, set<string*>& ignor_taxa) {
    if (leaf->left == 0 && leaf->right == 0) {
	if (taxa.find(leaf->nodelabel)!=taxa.end()) return 'T';
	else if (ignor_taxa.find(leaf->nodelabel)!=ignor_taxa.end()) return 'I';
	else return 'F';
    }
    char left = 'I';
    char right = 'I';
    if (leaf->left != 0) left = is_monophyletic_unrooted(leaf->left, taxa, ignor_taxa);
    if (leaf->right != 0) right = is_monophyletic_unrooted(leaf->right, taxa, ignor_taxa);
    
    if (left == 'X' || right == 'X') return 'X';
    else if (left == 'I' && right == 'I') return 'I';
    else if (left != 'I' && right == 'I') return left;
    else if (left == 'I' && right != 'I') return right;
    else if (left == 'T' && right == 'T') return 'T';
    else if (left == 'F' && right == 'F') return 'F';
    else if ((left == 'F' && right == 'T') || (left == 'T' && right == 'F')) return 'C';
    else if ((left == 'D' && right == 'T') ||  (left == 'T' && right == 'D')) return 'D';
    else if ((left == 'E' && right == 'F') ||  (left == 'F' && right == 'E')) return 'E';
    else if ((left == 'C' && right == 'T') ||  (left == 'T' && right == 'C')) return 'D';
    else if ((left == 'C' && right == 'F') ||  (left == 'F' && right == 'C')) return 'E';
    return 'X'; // is not monophyletic
}

char tree::is_monophyletic_unrooted (node* leaf, set<string*>& taxa) {
    if (leaf->left == 0 && leaf->right == 0) {
        if (taxa.find(leaf->nodelabel)!=taxa.end()) return 'T';
        else return 'F';
    }
    char left = 'T';
    char right = 'T';
    if (leaf->left != 0) left = is_monophyletic_unrooted(leaf->left, taxa);
    if (leaf->right != 0) right = is_monophyletic_unrooted(leaf->right, taxa);

    if (left == 'X' || right == 'X') return 'X';
    else if (left == 'T' && right == 'T') return 'T';
    else if (left == 'F' && right == 'F') return 'F';
    else if ((left == 'F' && right == 'T') || (left == 'T' && right == 'F')) return 'C';
    else if ((left == 'D' && right == 'T') ||  (left == 'T' && right == 'D')) return 'D';
    else if ((left == 'E' && right == 'F') ||  (left == 'F' && right == 'E')) return 'E';
    else if ((left == 'C' && right == 'T') ||  (left == 'T' && right == 'C')) return 'D';
    else if ((left == 'C' && right == 'F') ||  (left == 'F' && right == 'C')) return 'E';
    return 'X'; // is not monophyletic
}

tree::node* tree::max_proportion_taxa(node* leaf, const string& taxa) {
    if ( taxa.empty() ) return 0;
    set<string*> query_taxa;
    int unsigned length = taxa.length();
    string taxon;
    for (int unsigned i=0; i <= length; ++i) {
        if (taxa[i]==',' || i==length) {
            string* query = nodelabels.find_string(taxon);
            if (query != 0) query_taxa.insert(query);
            taxon.clear();
        }
        else taxon += taxa[i];
    }
    unsigned int in_taxa_above(0);
    unsigned int out_taxa_above(0);
    double max_proportion(0.0);
    unsigned int tot_n_taxta = n_sub_tips(leaf);
    #ifdef DEBUG
    cerr << "N taxa: " << query_taxa.size() /*<< " (" << taxa << ')'*/ << endl;
    #endif //DEBUG
    return max_proportion_taxa(leaf, query_taxa, in_taxa_above, out_taxa_above, tot_n_taxta, max_proportion);
}

tree::node* tree::max_proportion_taxa(node* leaf, const set<string*>& taxa, unsigned int& in_taxa_above, unsigned int& out_taxa_above, const unsigned int tot_n_taxa, double& max_proportion) {
    if (leaf == 0)  return 0;
    if (leaf->right == 0 && leaf->left == 0) {
	if (taxa.find(leaf->nodelabel) != taxa.end()) {
	    ++in_taxa_above;
	    if (1.0/taxa.size() > max_proportion) {
		max_proportion = 1.0/taxa.size();
		return leaf;
	    }
	}
	else {
	    ++out_taxa_above;
	}
	return 0;
    }
    else {
	unsigned int in_taxa_above_right(0);
	unsigned int out_taxa_above_right(0);
	node* left_return(0);
	node* right_return(0);
	if (leaf->left != 0) left_return = max_proportion_taxa(leaf->left, taxa, in_taxa_above, out_taxa_above, tot_n_taxa, max_proportion);
	if (leaf->right != 0) right_return = max_proportion_taxa(leaf->right, taxa, in_taxa_above_right, out_taxa_above_right, tot_n_taxa, max_proportion);
	in_taxa_above += in_taxa_above_right;
	out_taxa_above += out_taxa_above_right;
	if (in_taxa_above > out_taxa_above && (double)(in_taxa_above-out_taxa_above)/taxa.size() > max_proportion) {
	    max_proportion = (double)(in_taxa_above-out_taxa_above)/taxa.size();
	    #ifdef DEBUG
	    cerr << "In above: " << in_taxa_above << " Out above: " << out_taxa_above << " Prop above: " << (double)(in_taxa_above-out_taxa_above)/taxa.size() << endl;
	    cerr << "Max new up: " << max_proportion << endl;
	    #endif //DEBUG
	    return leaf;
	}
	else if ((taxa.size()-in_taxa_above) > (tot_n_taxa-taxa.size()-out_taxa_above) && (double)((taxa.size()-in_taxa_above-tot_n_taxa+taxa.size()+out_taxa_above))/taxa.size() > max_proportion) {
	    max_proportion = (double)((taxa.size()-in_taxa_above-tot_n_taxa+taxa.size()+out_taxa_above))/taxa.size();
	    #ifdef DEBUG
	    cerr << "In above: " << in_taxa_above << " Out above: " << out_taxa_above << " Prop below: " << (double)((taxa.size()-in_taxa_above-tot_n_taxa+taxa.size()+out_taxa_above))/taxa.size() << ' ' << (double)((2*taxa.size()-in_taxa_above-(tot_n_taxa-out_taxa_above))/taxa.size()) << endl;
	    cerr << "Max new down: " << max_proportion << endl;
	    #endif //DEBUG
	    return leaf;
	}
	else if (right_return != 0) return right_return;
	else return left_return;
    }
}

void tree::print_relaxed_outgroup( const string& taxa ) {
    node* theNode = max_proportion_taxa( root, taxa );
    if (theNode != 0) {
	if (theNode == root) {
	    #ifdef DEBUG
	    cerr << "The root is the node to root at" << endl;
	    #endif // DEBUG
	    if (n_sub_tips(theNode->left) > n_sub_tips(theNode->right)) print_tips(theNode->right,"",",");
	    else print_tips(theNode->left, "", ",");
	}
	else if (n_sub_tips(theNode) > n_sub_tips(root)/2) {
	    #ifdef DEBUG
	    cerr << "The non descendants are the smaller group." << endl;
	    #endif // DEBUG
	    print_non_descendants(cout, theNode->parent, theNode, "", ",");
	}
	else {
	    #ifdef DEBUG
	    cerr << "The descendants are the smaller group." << endl;
	    #endif // DEBUG
	    print_tips (theNode, "", ",");
	}
    }
}

bool tree::is_nested_in ( const node* ancestor, const node* descendent) {
    if (ancestor == 0 || descendent == 0) return false;
    if (ancestor == descendent) return true;
    bool answerLeft(false);
    bool answerRight(false);
    if (ancestor->left != 0) answerLeft = is_nested_in( ancestor->left, descendent);
    if (ancestor->right != 0) answerRight = is_nested_in( ancestor->right, descendent);
    return answerLeft || answerRight;
}

char tree::get_conflict_nodes (node* leaf, set<string*>& taxa, set<string*>& ignor_taxa, set<node*>& conflict_nodes) {
    if (leaf == 0) return 'I';
    if (leaf->left == 0 && leaf->right == 0) {
	if (ignor_taxa.find(leaf->nodelabel)!=ignor_taxa.end()) return 'I';
	else if (taxa.find(leaf->nodelabel)!=taxa.end()) return 'T';
	else return 'A';
    }
    char left = get_conflict_nodes(leaf->left, taxa, ignor_taxa,conflict_nodes);
    char right = get_conflict_nodes(leaf->right, taxa, ignor_taxa,conflict_nodes);
    if (left == 'I' && right == 'I') return 'I';
    else if (left != 'I' && right == 'I') return left;
    else if (left == 'I' && right != 'I') return right;
    else if (left == 'T' && right == 'T') return 'T';
    else if (left == 'A' && right == 'A') return 'A';
    else if ((left == 'A' && right == 'T') || (left == 'T' && right == 'A')) conflict_nodes.insert(leaf);
    return 'C';
}

void tree::tips_not_present_in_set (node* leaf, set<string*>& taxa, set<string*>& output) {
    if (leaf!=0) {
        if (leaf->left == 0 && leaf->right == 0) {
            if (taxa.find(leaf->nodelabel)==taxa.end()) {
                output.insert(leaf->nodelabel);
            }
        }
        else {
            if (leaf->left != 0) tips_not_present_in_set(leaf->left, taxa, output);
            if (leaf->right != 0) tips_not_present_in_set(leaf->right, taxa, output);
        }
    }
}

string tree::supported_not_monophyletic ( set<string*>& taxa, const float supp, set<string*>& ignor_taxa ) {
    set<string*> conf_taxa;
    if (taxa.size()>1) {
	node* ancestor = most_recent_common_ancestor( taxa );
	if (ancestor == 0) {
	    std::cerr << "ERROR!!! Could not find most recent common ancestor for ";
	    for (set<string*>::iterator i=taxa.begin(); i!=taxa.end();++i) {
		if (i != taxa.begin()) std::cerr << ", ";
		std::cerr << *(*i);
	    }
	    std::cerr << '.' << std::endl;
	}
	else {
	    set<node*> conflict_nodes;
	    get_conflict_nodes(ancestor->left, taxa, ignor_taxa, conflict_nodes); // find conflicting taxa to the left
	    get_conflict_nodes(ancestor->right, taxa, ignor_taxa, conflict_nodes); // find conflicting taxa to the right
	    if (!conflict_nodes.empty() && !(ancestor==root && conflict_nodes.size()<2)) {
		set<node*> supported_inlucive_clades;
	    	for (set<node*>::const_iterator i=conflict_nodes.begin(); i!=conflict_nodes.end(); ++i) {
    		    node* present=*i;
		    node* basal_supported=0;
		    while (present!=0 && present != ancestor) {
			if (supported_inlucive_clades.find(present) != supported_inlucive_clades.end()) break;
			if (present->nodelabel != 0 && strtod(present->nodelabel->c_str(),0) > supp) basal_supported=present;
			present = present->parent;
		    }
		    if (basal_supported != 0) supported_inlucive_clades.insert(basal_supported);
		}
		set<string*> no_conflict_taxa = taxa;
		for (set<string*>::iterator i=ignor_taxa.begin(); i != ignor_taxa.end(); ++i) {
		    no_conflict_taxa.insert(*i);
		}
		for (set<node*>::iterator i=supported_inlucive_clades.begin(); i != supported_inlucive_clades.end(); ++i) {
		    tips_not_present_in_set (*i, no_conflict_taxa, conf_taxa);
		}
	    }
	}
    }
    string return_string;
    for (set<string*>::iterator i=conf_taxa.begin(); i!=conf_taxa.end();++i) {
	if (i!=conf_taxa.begin()) return_string+=',';
	return_string+=*(*i);
    }
    return return_string;
}

// Print clades where there is a supported conflict between trees, possibly as SVG
void tree::print_conflict_clades ( node* leaf, tree* tree, const float supp, set<string*>& not_in_tree, set<string*>& not_in_comp_tree, bool svg ) {
    if (leaf != 0) {
        if ( leaf->left != 0 || leaf->right != 0) {
            if (leaf->left != 0) print_conflict_clades (leaf->left, tree, supp, not_in_tree, not_in_comp_tree, svg);
            if (leaf->right != 0) print_conflict_clades (leaf->right, tree, supp, not_in_tree, not_in_comp_tree, svg);
	    if (leaf->nodelabel != 0) {
		float node_support = strtod(leaf->nodelabel->c_str(),0);
		if (node_support > supp) {
		    set<string*> tip_names;
		    this->tips_not_present_in_set( leaf, not_in_comp_tree, tip_names ); // get tips that are in the clade and in tree that is compared too
		    string missfits = tree->supported_not_monophyletic(tip_names, supp, not_in_tree);
		    if (!missfits.empty()) {
			string taxa_names;
			for (set<string*>::iterator i=tip_names.begin(); i!=tip_names.end(); ++i) {
			    if (i != tip_names.begin()) taxa_names+=',';
			    taxa_names+=*(*i);
			}
			if (svg) std::cout << "<p>" << endl;
			cout << "The taxa: " << missfits << ", are supported as nested in " << taxa_names << endl;
			if (svg) {
			    std::cout << "</p>" << endl;
			    std::cout << "<h2>First tree:</h2>" << endl;
			    string tip_color = "red(" + missfits + ")green(" + taxa_names + ")";
			    print_svg_autoscale_no_head(false,true,10.0,tip_color, true);
			    std::cout << "<h2>Second tree:</h2>" << endl;
			    tree->print_svg_autoscale_no_head(false,true,10.0,tip_color, true);
			}
		    }
		}
	    }
        }
    }
}
// Function to find taxa that is not present in other tree and to call function to print conflicts between trees
void tree::print_conflict_clades_reduced ( tree* tree, const float supp, bool svg ) {
    set<string*> taxa_tree1; tip_names( root, taxa_tree1 );
    set<string*> taxa_tree2; tree->tip_names( tree->root, taxa_tree2 );
    set<string*> exclude_taxa1;
    for (set<string*>::iterator i=taxa_tree1.begin(); i != taxa_tree1.end(); ++i) {
	if (taxa_tree2.find(*i)==taxa_tree2.end()) exclude_taxa1.insert(*i);
    }
    set<string*> exclude_taxa2;
    for (set<string*>::iterator i=taxa_tree2.begin(); i != taxa_tree2.end(); ++i) {
	if (taxa_tree1.find(*i)==taxa_tree1.end()) exclude_taxa2.insert(*i);
    }
    print_conflict_clades ( root, tree, supp, exclude_taxa2, exclude_taxa1, svg );
}

unsigned int tree::splits_not_in_B(node* leaf, tree* B, set<string*>& split, const unsigned int n_shared_taxa, set<string*>& not_in_A, set<string*>& not_in_B) {
    if (leaf == 0 || (leaf->left == 0 && leaf->right==0)) {
	if (not_in_B.find(leaf->nodelabel)==not_in_B.end()) split.insert(leaf->nodelabel);
	return 0;
    }
    unsigned int score(0);
    set<string*> split_left;
    set<string*> split_right;
    if (leaf->left !=0) score += splits_not_in_B(leaf->left, B, split_left, n_shared_taxa, not_in_A, not_in_B);
    if (leaf->right !=0) score += splits_not_in_B(leaf->right, B, split_right, n_shared_taxa, not_in_A, not_in_B);
    if (!split_left.empty() && !split_right.empty()) {
	unsigned int left_size = split_left.size();
	if (left_size>1) {
	    char mono = B->is_monophyletic_unrooted(B->root,split_left);//,not_in_A);
	    if (mono == 'X') ++score;
	}
	unsigned int right_size = split_right.size();
	if ( right_size>1 && left_size+right_size < n_shared_taxa ) {
	    char mono = B->is_monophyletic_unrooted(B->root,split_right);//,not_in_A);
	    if (mono == 'X') ++score;
	}
    }
    split.insert(split_left.begin(),split_left.end());
    split.insert(split_right.begin(),split_right.end());
    return score;
}

unsigned int tree::robinson_fould(tree* B) {
    set<string*> taxa_tree1; tip_names( root, taxa_tree1 );
    set<string*> taxa_tree2; B->tip_names( B->root, taxa_tree2 );
    set<string*> exclude_taxa1;
    for (set<string*>::iterator i=taxa_tree1.begin(); i != taxa_tree1.end(); ++i) {
        if (taxa_tree2.find(*i)==taxa_tree2.end()) exclude_taxa1.insert(*i);
    }
    set<string*> exclude_taxa2;
    for (set<string*>::iterator i=taxa_tree2.begin(); i != taxa_tree2.end(); ++i) {
        if (taxa_tree1.find(*i)==taxa_tree1.end()) exclude_taxa2.insert(*i);
    }
    set<string*> complete_split;
    unsigned int n_shared_taxa = shared_tips(B);
    unsigned int return_value(0);
    return_value += splits_not_in_B(root,B, complete_split, n_shared_taxa, exclude_taxa2, exclude_taxa1);
    complete_split.clear();
    return_value += B->splits_not_in_B(B->root, this, complete_split, n_shared_taxa, exclude_taxa1, exclude_taxa2);

    return return_value;
}

void tree::add_to_support(node* leaf, tree* B, set<string*>& split, const bool rooted) {
    if (leaf!=0) {
    	if (leaf->left==0 && leaf->right==0) {
    	    if (leaf->nodelabel!=0) split.insert(leaf->nodelabel);
    	}
    	else {
	    set<string*> descendants;
	    if (leaf->left!=0) add_to_support(leaf->left, B, descendants, rooted);
	    if (leaf->right!=0) add_to_support(leaf->right, B, descendants, rooted);
	    if (leaf != root) {
		char mono('X');
		if (rooted) {
		   if (B->is_monophyletic(B->root,descendants)) mono='T';
		}
	       	else {
		    if (!(leaf->parent == root && n_sub_tips(root)-n_sub_tips(leaf) == 1)) mono = B->is_monophyletic_unrooted(B->root,descendants);
		    #ifdef DEBUG
		    cerr << "N descendants: " << descendants.size() << endl;
		    cerr << "Mono code: " << mono << endl;
		    #endif //DEBUG
		}
		if (mono != 'X' && mono != 'F') {
		    stringstream converter;
		    int present_support(0);
		    if (leaf->nodelabel!=0) {
			present_support = atoi(leaf->nodelabel->c_str());
		    }
		    ++present_support;
		    converter << present_support;
		    string value(converter.str());
		    leaf->nodelabel = nodelabels.add_string(value);
		}
		split.insert(descendants.begin(),descendants.end());
	    }
      	}
    }
}

void tree::add_to_support(tree* B, const bool rooted) {
    set<string*> taxa_in_tree;
    add_to_support(root,B,taxa_in_tree, rooted);
}

void tree::add_one_to_each_support(node* leaf, const bool rooted) {
    if (leaf != 0) {
	if (leaf->left != 0) add_one_to_each_support(leaf->left, rooted);
	if (leaf->right != 0) add_one_to_each_support(leaf->right, rooted);
	if (leaf->left != 0 && leaf->right != 0 && leaf != root) {
	    if (rooted || !(leaf->parent == root && n_sub_tips(root)-n_sub_tips(leaf) == 1)) {
		stringstream converter;
		int present_support(0);
		if (leaf->nodelabel!=0) {
		    present_support = atoi(leaf->nodelabel->c_str());
		}
		++present_support;
		converter << present_support;
		string value(converter.str());
		leaf->nodelabel = nodelabels.add_string(value);
	    }
	}
    }
}

unsigned int tree::tips_present_in_set (node* leaf, set<string*>& taxa) {
    if (leaf==0) return 0;
    if (leaf->left==0 && leaf->right==0) {
	if (taxa.find(leaf->nodelabel)!=taxa.end()) return 1;
    }
    unsigned int left=0;
    unsigned int right=0;
    if (leaf->left != 0) left = tips_present_in_set(leaf->left, taxa);
    if (leaf->right != 0) right = tips_present_in_set(leaf->right, taxa);
    return left+right;
}

void tree::tips_not_shared( tree* B, const string& name_tree1, const string & name_tree2 ) {
    set<string*> taxa_tree1; tip_names( root, taxa_tree1 );
    set<string*> taxa_tree2; B->tip_names( B->root, taxa_tree2 );
    bool printed = false;
    for (set<string*>::iterator i=taxa_tree1.begin(); i != taxa_tree1.end(); ++i) {
        if (taxa_tree2.find(*i)==taxa_tree2.end()) {
	    if (printed) std::cout << ", ";
	    std::cout << *(*i); printed = true;
	}
    }
    if (printed) std::cout << " in " << name_tree1 << " are missing from " << name_tree2 << endl;
    else std::cout << "All tips in "  << name_tree1 << " are present in " << name_tree2 << endl;
}

unsigned int tree::shared_tips(tree* B) {
    set<string*> taxa; tip_names( root, taxa );
    return B->tips_present_in_set(B->root, taxa);
}

void tree::adjustedMPL (const string nodes_and_ages, const int n_char){
    map<node*,double> given_nodeages;
    string taxa;
    string age;
    char mode='T';
    int length=nodes_and_ages.length();
    for (int i=0; i<=length;++i) {
	if (i==length || nodes_and_ages[i]==';') {
	    if (!age.empty() && !taxa.empty()) {
		double given_age=atof(age.c_str());
		if (given_age<=0) std::cerr << "Suspicious age for node " << taxa << ": " << given_age << std::endl;
		node* MRCA=0;
		if (taxa.compare("root") ==0 ) MRCA = root;
	       	else MRCA = most_recent_common_ancestor(taxa);
		if (MRCA != 0) {
		    given_nodeages[MRCA] = given_age;
		}
	    }
	    else std::cerr << "Error in parsing age of clades." << std::endl; 
	    taxa.clear();
	    age.clear();
	    mode='T';
	}
	else if (nodes_and_ages[i]==':') mode='A';
	else if (mode == 'T') taxa+=nodes_and_ages[i];
	else if (mode == 'A') age+=nodes_and_ages[i];
    }
    adjustedMPL(given_nodeages,n_char);
}

void tree::adjustedMPL(map<node*,double>& given_nodeages, unsigned int n_char) {
    vector<double> rates;
    for (map<node*,double>::iterator i=given_nodeages.begin(); i != given_nodeages.end(); ++i) {
	rates.push_back(local_adjustedMPL(i->first, given_nodeages, n_char));
	// Add rate to nodelabel for nodes with ages given
	std::stringstream converter;
	converter << "[&rate=" << rates.back() << ']';
	string temp;
	if (i->first->nodelabel != 0) temp=*(i->first->nodelabel);
	temp += converter.str();
 	i->first->nodelabel = nodelabels.add_string(temp);
    }
    if (given_nodeages.find(root) == given_nodeages.end()) {
	if (!rates.empty()) {
	    sort(rates.begin(),rates.end());
	    double root_rate=rates[rates.size()/2];
	    int_double2 values_left = calculate_mean_path_root(root->left,given_nodeages);
	    int_double2 values_right = calculate_mean_path_root(root->right,given_nodeages);
	    double mean_path = (values_left.c+values_right.c+root_rate*n_char*(values_left.a+values_right.a)) / (values_left.n+values_right.n);
	    double root_age = mean_path/(n_char*root_rate);
	    given_nodeages[root]=root_age;
	    map<node*,double> node_depths;
	    if (root->left!=0) calculate_node_depth(root->left,root_age,mean_path,root_rate,given_nodeages,node_depths);
	    if (root->right!=0) calculate_node_depth(root->right,root_age,mean_path,root_rate,given_nodeages,node_depths);
	    assign_branches_based_on_nodeheights(root,node_depths,given_nodeages);
	    std::stringstream converter;
	    converter << "[&rate=" << root_rate << ']';
	    string temp;
	    if (root->nodelabel != 0) temp=*(root->nodelabel);
	    temp += converter.str();
	    root->nodelabel = nodelabels.add_string(temp);
	}
	else std::cerr << "Warning!!! No rates available to estimate root age. Branches with the root age are not rescaled!!!" << std::endl;
    }
}

tree::int_double2 tree::calculate_mean_path_root( node* leaf, map<node*,double>& given_nodeages ) {
    int_double2 values_left = {0,0.0,0.0};
    if (leaf==0) return values_left;
    else {
	if (given_nodeages.find(leaf) != given_nodeages.end()) {
	    values_left.n = n_sub_tips(leaf);
	    values_left.a = values_left.n*given_nodeages[leaf];
	    values_left.c = leaf->branchlength*values_left.n;
	    return values_left;
	}
	if (leaf->left==0 && leaf->right==0) {
	    values_left.n = 1;
	    values_left.c = leaf->branchlength;
            return values_left;
	}
	if (leaf->left!=0) values_left = calculate_mean_path_root(leaf->left, given_nodeages);
	int_double2 values_right = {0,0.0,0.0};
	if (leaf->right!=0) values_right = calculate_mean_path_root(leaf->right, given_nodeages);
	values_left.n += values_right.n;
	values_left.a += values_right.a;
	values_left.c += values_right.c+values_left.n*leaf->branchlength;
	return values_left;
    }
}

double tree::local_adjustedMPL(node* leaf, map<node*,double>& given_nodeages,unsigned int n_char) {
    int_double2 values_left = {0,0.0,0.0};
    if (leaf->left!=0) values_left = calculate_mean_path(leaf->left,given_nodeages[leaf],given_nodeages);
    int_double2 values_right = {0,0.0,0.0};
    if (leaf->right!=0) values_right = calculate_mean_path(leaf->right,given_nodeages[leaf],given_nodeages);
    double mean_path = (values_left.c+values_right.c)/n_sub_tips(leaf);
    double rate = mean_path/(n_char*given_nodeages[leaf]);
    map<node*,double> node_depths;
    if (leaf->left!=0) calculate_node_depth(leaf->left,given_nodeages[leaf],mean_path,rate,given_nodeages,node_depths);
    if (leaf->right!=0) calculate_node_depth(leaf->right,given_nodeages[leaf],mean_path,rate,given_nodeages,node_depths);
    assign_branches_based_on_nodeheights(leaf,node_depths,given_nodeages);
    return rate;
}

tree::int_double2 tree::calculate_mean_path( node* leaf, const double root_age, map<node*,double>& given_nodeages) {
    int_double2 values_left={0,0.0,0.0};
    if (leaf==0) return values_left;
    else {
	if (given_nodeages.find(leaf) != given_nodeages.end()) {
	    int n=n_sub_tips(leaf);
	    values_left.a=n*(root_age/(root_age-given_nodeages[leaf]));
	    values_left.c=values_left.a;
	    return values_left;
	}
	else if (leaf->left==0 && leaf->right==0) {
	    values_left.n = 1;
	    values_left.c = leaf->branchlength;
	    return values_left;
	}
	if (leaf->left!=0) values_left=calculate_mean_path(leaf->left,root_age,given_nodeages);
	int_double2 values_right={0,0.0,0.0};
	if (leaf->right!=0) values_right=calculate_mean_path(leaf->right,root_age,given_nodeages);
	values_left.n += values_right.n;
	values_left.a += values_right.a;
	values_left.c += values_right.c+(values_left.n+values_left.a)*leaf->branchlength;//(values_left.n+values+values_left.a+values_right+a;
	return values_left;
    }
}

void tree::assign_branches_based_on_nodeheights(node* leaf, map<node*,double>& node_depths, map<node*,double> given_nodeages) {
    for (map<node*,double>::iterator i=node_depths.begin(); i != node_depths.end(); ++i) {
        if (i->first->parent == leaf) {
	    i->first->branchlength = given_nodeages[leaf] - i->second;
	    #ifdef DEBUG
	    cerr << "A " << given_nodeages[leaf] << " " << i->second << " " << i->first->branchlength << endl;
	    #endif //DEBUG
	}
        else {
	    i->first->branchlength = node_depths[i->first->parent] - i->second;
	    #ifdef DEBUG
	    cerr << "B " << node_depths[i->first->parent] << " " << i->second << " " << i->first->branchlength << endl;
	    #endif //DEBUG
	}
    }
}

tree::int_double2 tree::calculate_node_depth(node* leaf, const double root_age, const double path_to_root, const double rate, map<node*,double>& given_nodeages, map<node*,double>& calculated_node_depths) {
    int_double2 values_left={0,0.0,0.0};
    if (leaf==0) return values_left;
    else {
	if (given_nodeages.find(leaf)!=given_nodeages.end()) {
	    values_left.n=n_sub_tips(leaf);
	    values_left.a=given_nodeages[leaf]*values_left.n;
	    values_left.c=leaf->branchlength*values_left.n;
	    calculated_node_depths[leaf]=given_nodeages[leaf];
	    return values_left;
	}
	if (leaf->left!=0) {
	    values_left=calculate_node_depth(leaf->left, root_age, path_to_root, rate, given_nodeages, calculated_node_depths);
	}
	int_double2 values_right={0,0.0,0.0};
	if (leaf->right!=0) {
	    values_right=calculate_node_depth(leaf->right, root_age, path_to_root, rate, given_nodeages, calculated_node_depths);
	}
    	if (leaf->left==0 && leaf->right==0) values_left.n=1;
	else values_left.n += values_right.n;
	double path = ((values_left.c+values_right.c)+rate*(values_left.a+values_right.a))/values_left.n;
	values_left.c += values_right.c+leaf->branchlength*values_left.n;
	values_left.a += values_right.a;
	calculated_node_depths[leaf]=root_age*(path/path_to_root);
	return values_left; 
    }
}

tree::int_double2 tree::test_clock_likness(node* leaf, const int n_char) {
    int_double2 values_left={0,0.0,0.0};
    if (leaf!=0) {
	if (leaf->left!=0) values_left=test_clock_likness(leaf->left, n_char);
	int_double2 values_right={0,0.0,0.0};
	if (leaf->right!=0) values_right=test_clock_likness(leaf->right, n_char);
	if (leaf->left==0 && leaf->right==0) values_left.n=1;	
	else {
	    double z = ((values_left.a/values_left.n) - (values_right.a/values_right.n))/sqrt((values_left.c/(values_left.n*values_left.n)) + (values_right.c/(values_right.n*values_right.n)));
	    std::stringstream convert;
	    convert << "[&clocklikness=" << z << ']';
	    leaf->nodelabel = nodelabels.add_string(convert.str());
	    values_left.n+=values_right.n;
	}
	values_left.a += values_right.a+(values_left.n*leaf->branchlength*n_char);
	values_left.c += values_right.c+(values_left.n*values_left.n*leaf->branchlength*n_char);
    }
    return values_left;
}

void tree::internal_nodes_stat( const string& type, bool sum, bool prod, bool average, bool median, bool SD, bool n_values_above, double cut_off, bool max, bool min, bool N, bool include_root ) {
    vector<double> values;
    if (include_root && root->nodelabel!=0) values.push_back(atof(root->nodelabel->c_str()));
    if (root->left!=0 && (root->left->left != 0 || root->left->right != 0)) get_internal_node_values(root->left, type, values);
    if (root->right!=0 && (root->right->left != 0 || root->right->right != 0)) get_internal_node_values(root->right, type, values);
    if (!values.empty()) {
	double sum=0;
	double prod =1;
	sort(values.begin(),values.end());
	unsigned int n_above=0;
	for (vector<double>::iterator i=values.begin(); i != values.end(); ++i) {
	    sum += *i;
	    prod *= *i;
	    if (*i > cut_off) ++n_above;
	}
	if (sum) std::cout << "Sum: " << sum << endl;
	if (prod) std::cout << "Prod: " << prod << endl;
	if (average) std::cout << "Average: " << sum/values.size() << endl;
	if (median) {
	    if (values.size()%2) std::cout << "Median: " << values[values.size()/2] << endl;
	    else std::cout << "Median: " << (values[values.size()/2]+values[(values.size()/2)+1])/2 << endl;
	}
	if (SD) {
	    double Var = 0;
	    double average = sum/values.size();
	    for (vector<double>::iterator i=values.begin(); i != values.end(); ++i) {
		Var += ((*i)-average) * ((*i)-average);
	    }
	    Var /= values.size();
	    std::cout << "Standard deviation: " << sqrt(Var) << endl;
	}
	if (n_values_above) std::cout << n_above << " values were above " << cut_off << endl;
	if (max) std::cout << "Max value: " << values[values.size()-1] << endl;
	if (min) std::cout << "Min value: " << values[0] << endl;
	if (N) std::cout << "N internal nodes with labels: " << values.size() << endl;
    }
    else if (type.empty()) std::cerr << "Found no node labels." << endl;
    else std::cerr << "Found no values for " << type << '.' << std::endl;
}

void tree::get_internal_node_values(node* leaf, const string& type, vector<double>& values) {
    if (leaf->left!=0 && (leaf->left->left!=0 || leaf->left->right!=0)) get_internal_node_values(leaf->left, type, values);
    if (leaf->right!=0 && (leaf->right->left!=0 || leaf->right->right!=0)) get_internal_node_values(leaf->right, type, values);
    if (leaf->nodelabel!=0) 
	add_value_from_label_to_vector (*leaf->nodelabel, type, values);
}

void tree::add_value_from_label_to_vector (const string& label, const string& param, vector<double>& values) {
    if (param.empty()) values.push_back(atof(label.c_str()));
    else {
	std::size_t found = label.find("[&");
	found = label.find(param,found);
	while (found < label.length() && label[found] != '=') ++found;
	if (found < label.length()) {
	    ++found;
	    std::size_t end=found;
	    while (end < label.length() && label[end] != ',' && label[end] != ']') ++end;
	    if (end == label.length()) end = found;
	    //else --end;
	    if(end-found>0) {
		values.push_back(atof(label.substr(found,end-found).c_str()));
		#ifdef DEBUG
		cerr << "Lable value: " << label.substr(found,end-found) << " " << atof(label.substr(found,end-found).c_str()) << endl;
		#endif //DEBUG
	    }
	}
    }
}

void tree::code_clades_in_matrix(node* leaf, map<string,vector<char> >& matrix, set<string>& absent_taxa) {
    if (leaf== 0 || ( leaf->left == 0 && leaf->right == 0 )) return;
    if (leaf->left != 0) code_clades_in_matrix(leaf->left, matrix, absent_taxa);
    if (leaf->right != 0) code_clades_in_matrix(leaf->right, matrix, absent_taxa);
    set<string*> one_side;
    tip_names(leaf, one_side);
    if (one_side.size() > 1) {
	set<string*> other_side;
	tips_not_present_in_set (root, one_side, other_side);
	if (other_side.size() > 1) {
	    for (set<string*>::const_iterator i=one_side.begin(); i!= one_side.end(); ++i) {
		matrix[*(*i)].push_back('1');
	    }
	    for (set<string*>::const_iterator i=other_side.begin(); i!= other_side.end(); ++i) {
		matrix[*(*i)].push_back('0');
	    }
	    for (set<string>::const_iterator i=absent_taxa.begin(); i!= absent_taxa.end(); ++i) {
		matrix[*i].push_back('-');
	    }
	}
    }
}

void tree::add_taxa_to_matrix(node* leaf, map<string,vector<char> >& matrix, unsigned int length) {
    if (leaf==0) return;
    if (leaf->left==0 && leaf->right==0) {
	map<string,vector<char> >::iterator match = matrix.find(*leaf->nodelabel);
	if (match == matrix.end()) {
	    matrix[*leaf->nodelabel] = vector<char>();
	    match = matrix.find(*leaf->nodelabel);
 	    for (unsigned int n=0; n<length; ++n) match->second.push_back('-');
	}
    }
    else {
	if (leaf->left!=0) add_taxa_to_matrix(leaf->left,matrix,length);
	if (leaf->right!=0) add_taxa_to_matrix(leaf->right,matrix,length);
    }
}

void tree::add_to_matrix_representation(map<string,vector<char> >& matrix) {
    set<string> absent_taxa;
    for (map<string,vector<char> >::const_iterator i=matrix.begin();i!=matrix.end();++i) {
	if(find_taxon_tip(root, i->first)==0) absent_taxa.insert(i->first);
    }
    if (!matrix.empty()) {
	unsigned int length_matrix = matrix.begin()->second.size();
	add_taxa_to_matrix(root,matrix,length_matrix);
    }
    if (root->left != 0) code_clades_in_matrix(root->left, matrix, absent_taxa);
    if (root->right != 0) {
	if (root->right->left != 0) code_clades_in_matrix(root->right->left, matrix, absent_taxa);
	if (root->right->right != 0) code_clades_in_matrix(root->right->right, matrix, absent_taxa);
    }
}

void tree::copy (node* B, const node* A, node* parent) {
    B->parent = new node;
    B->parent = parent;
    B->branchlength = A->branchlength;
    B->nodelabel = A->nodelabel;
    if (A->left != 0) {
	B->left = new node;
	copy(B->left, A->left, B);
    }
    else B->left=0;
    if (A->right != 0) {
	B->right = new node;
	copy(B->right,A->right, B);
    }
    else B->right = 0;
}
unsigned int tree::get_internal_node_by_number(node* leaf, node*& tip, unsigned int node_no) {
    if (leaf == 0 || (leaf->left==0 && leaf->right==0)) return node_no;
    --node_no;
    if (node_no == 0) { tip = leaf; return node_no; }
    else {
	if (leaf->left != 0) node_no = get_internal_node_by_number(leaf->left, tip, node_no);
	if (tip == 0 && leaf->right !=0) node_no = get_internal_node_by_number(leaf->right, tip, node_no);
	return node_no;
    }
}

tree::node* tree::get_internal_branch_by_number_unrooted(unsigned int branch_no) {
    node* return_node = 0;
    if (branch_no == 0) return return_node;
    else if (n_sub_tips(root->left) < 2)
	get_internal_node_by_number(root->right, return_node, branch_no);
    else if (n_sub_tips(root->right) < 2)
	get_internal_node_by_number(root->left, return_node, branch_no);
    else {
	++branch_no;
	get_internal_node_by_number(root, return_node, branch_no);
    }
    return return_node;
}

void tree::interchange(node* leaf, bool left) {
    if (leaf == 0 || leaf->parent == 0 || (leaf->left ==0 && leaf->right ==0)) return;
    node* swap_partner = 0; 
    if (leaf->parent->left != 0 && leaf->parent->left != leaf) {
	if (leaf->parent != root) 
	    swap_partner = leaf->parent->left;
	else if (leaf->parent->left->left != 0)
	    swap_partner = leaf->parent->left->left;
    }
    else if (leaf->parent->right != 0 && leaf->parent->right != leaf) {
	if (leaf->parent != root) 
	    swap_partner = leaf->parent->right;
	else if (leaf->parent->right->left != 0)
	    swap_partner = leaf->parent->right->left;
    }
    if (swap_partner != 0) {
	node* temp = 0;
	if (left && leaf->left !=0) {
	    temp = leaf->left;
	    leaf->left = swap_partner;
	}
	if (!left && leaf->right !=0) {
	    temp = leaf->right;
	    leaf->right = swap_partner;
	}
	if (temp != 0) {
	    temp->parent = swap_partner->parent;
	    if (temp->parent->left == swap_partner) {
		temp->parent->left = temp;
	    }
	    else if (temp->parent->right == swap_partner) {
                swap_partner->parent->right = temp;
            }
	    swap_partner->parent = leaf;
	}
    }
}

void tree::nni (unsigned int branch_no, tree& L, tree& R) {
    L = *this;
    node* swap_branch = 0;
    if (branch_no == 0) branch_no = (rand() % n_tips()) +1;
    swap_branch = L.get_internal_branch_by_number_unrooted(branch_no);
    #ifdef DEBUG
    cout << "N tips: " << L.n_sub_tips(root) -  L.n_sub_tips(swap_branch) << endl;
    #endif //DEBUG
    if ( branch_no == 1 && (L.n_sub_tips(swap_branch)==1  || L.n_sub_tips(L.root)-L.n_sub_tips(swap_branch)==1 )) {
	L.destroy_tree(L.root->left);
	L.root->left=0;
	L.destroy_tree(L.root->right);
	L.root->right=0;
	R.destroy_tree(R.root->left);
	R.root->left=0;
	R.destroy_tree(R.root->right);
	R.root->right=0;
	return;
    }
    if (swap_branch != 0) {
	L.interchange(swap_branch, true);
    }
    R = *this;
    swap_branch = 0;
    swap_branch = R.get_internal_branch_by_number_unrooted(branch_no);
    if (swap_branch != 0) {
	R.interchange(swap_branch, false);
    }
}

unsigned int tree::fitch_parsimony (node* leaf, map<node*, parsimony_character_vector>& characters, const unsigned int start_char, const unsigned int end_char) {
    if (leaf == 0) return 0;
    if (leaf->left == 0 && leaf->right == 0) {
	if (characters.find(leaf) == characters.end()) {
	    bitset<SIZE> character;
	    character.flip();
	    for (unsigned int i=start_char; i <= end_char; ++i) {
		characters[leaf].add_character(i,character); // if no state given for the tip, set as uncertain
	    }
	    #ifdef DEBUGTREEATOR
	    cerr << "No char for " << *leaf->nodelabel << endl;
	    #endif //DEBUGTREEATOR
	}
	return 0;
    }
    unsigned int score = 0;
    if (leaf->left != 0) score += fitch_parsimony(leaf->left, characters, start_char, end_char);
    if (leaf->right != 0) score += fitch_parsimony(leaf->right, characters, start_char, end_char);
    if (characters.find(leaf) == characters.end())
	characters[leaf] = parsimony_character_vector(); // if nothing there initiate it
    for (unsigned int i=start_char; i <= end_char; ++i) {
	//bitset<24>* this_node = characters[leaf].get_character(i);
	if (!characters[leaf].get_character(i).any()) { // if no state given calc state from descendants
	    if (leaf->left != 0 && characters.find(leaf->left) != characters.end() && leaf->right != 0 && characters.find(leaf->right) != characters.end()) { // if both descendentas are initiated
		bitset<SIZE> intersect = characters[leaf->left].get_character(i); 
		intersect &= characters[leaf->right].get_character(i);
		if (intersect.any()) characters[leaf].add_character(i,intersect); // set to intersect
		else {
		    intersect = characters[leaf->left].get_character(i); 
		    intersect |= characters[leaf->right].get_character(i);
		    characters[leaf].add_character(i,intersect); // if no intersect set to union
		    ++score; // add to score
		}
	    }
	    else if (leaf->left != 0 && characters.find(leaf->left) != characters.end()) characters[leaf].add_character(i,characters[leaf->left].get_character(i));
	    else if (leaf->right != 0 && characters.find(leaf->right) != characters.end()) characters[leaf].add_character(i,characters[leaf->right].get_character(i));
	    else { characters[leaf].get_character(i).flip();
		#ifdef DEBUG
		cerr << "ERROR!!! Descendants without characters" << endl;
		#endif //DEBUG
	    }
	}
	else { // if state already given
	    #ifdef DEBUG
    	    cerr << "ERROR!!! Character given for internal node" << endl;
	    #endif //DEBUG
	    if (leaf->left != 0 && characters.find(leaf->left) != characters.end()) {
		bitset<SIZE> intersect = characters[leaf].get_character(i);
		intersect &= characters[leaf->left].get_character(i);
		if (!intersect.any()) ++score; // if no intersect increase score
	    }
	    if (leaf->right != 0 && characters.find(leaf->right) != characters.end()) {
		bitset<SIZE> intersect = characters[leaf].get_character(i);
		intersect &= characters[leaf->right].get_character(i);
		if (!intersect.any()) ++score; // if no intersect increase score
	    }
	}
    }
    #ifdef DEBUG
    cerr << score << endl;
    #endif //DEBUG
    characters[leaf].set_score(score);
    return score;
}

unsigned int tree::fitch_parsimony (vector<character_vector>& characters, bool ancestral_states, bool get_branch_lengths, map<char, bitset<SIZE> >& alphabet) {
    map<node*, parsimony_character_vector > decoded_characters; // map characters to nodes
    unsigned int max_n_char = 0;
    for (vector<character_vector>::iterator i= characters.begin(); i != characters.end(); ++i) {
	node* taxa = find_taxon_tip(root, i->get_taxon());
	if (taxa == 0) std::cerr << i->get_taxon() << " is not present in tree." << std::endl;
	else {
	    decoded_characters[taxa] = parsimony_character_vector(*i);
	    unsigned int n_char = i->n_char();
	    if (max_n_char != 0 && max_n_char != n_char) std::cerr << "WARNING!!! " << *(taxa->nodelabel) << " has " << n_char << " chatacters while previous taxa has had up to " << max_n_char << " characters." << std::endl;
	    if (n_char > max_n_char) max_n_char = n_char;
	}
    }
    #ifdef DEBUG
    cerr << "N char: " << max_n_char << endl;
    #endif //DEBUG
    if (max_n_char==0) return 0;
    unsigned int score = fitch_parsimony(root, decoded_characters, 0, max_n_char-1);
    if (ancestral_states || get_branch_lengths) {
	unsigned int n = decoded_characters[root].n_char();
	parsimony_character_vector prefered; // prefered is used to select along what branch a character change
	bitset<SIZE> temp;
	for (unsigned int i=0; i<n; ++i) { // pick random possible root state for each character
	    temp = support_functions::pick_a_random_true_bit(decoded_characters[root].get_character(i));
	    prefered.add_character(i,temp);
	}
	fitch_parsimony_second_pass(root, decoded_characters, prefered, get_branch_lengths, ancestral_states, false, 0, max_n_char-1, alphabet);
    }
    return score;
}

void tree::fitch_parsimony_second_pass (node* leaf, map<node*, parsimony_character_vector > characters, parsimony_character_vector prefered, bool calc_branch_length, bool draw_ancestral_state, bool reconstruct_tip_state, const unsigned int start_char, const unsigned int end_char, map<char, bitset<SIZE> >& alphabet ) {
    if (leaf==0) return;
    if (leaf->parent !=0 && characters.find(leaf->parent) != characters.end()) {
	unsigned int branch(0);
	for (unsigned int i=start_char; i<=end_char; ++i) {
	    bitset<SIZE> preliminary = characters[leaf].get_character(i);
	    bitset<SIZE> temp = characters[leaf->parent].get_character(i);
	    temp ^= preliminary;
	    temp &= characters[leaf->parent].get_character(i);
	    // get final node set
	    if (temp.none() && (reconstruct_tip_state || leaf->left != 0 || leaf-> right != 0)) { // Fitch second round first option
		preliminary &= characters[leaf->parent].get_character(i);
		characters[leaf].add_character(i,preliminary);
	    }
	    else if (leaf->left != 0) {
		temp = characters[leaf->left].get_character(i);
		temp ^= preliminary; // get bits unique for any set
		temp &= preliminary; // get bits unique for preliminary
		//test if formed by union
		if (temp.any()) { // if any bits is uniqe for preliminary it was a union
		    preliminary |= characters[leaf->parent].get_character(i);
		    characters[leaf].add_character(i,preliminary); // add any states in ancestor
		}
		else if (leaf->right != 0) {
		    temp = characters[leaf->left].get_character(i) | characters[leaf->right].get_character(i);
		    temp &= characters[leaf->parent].get_character(i);
		    preliminary |= temp;
		    characters[leaf].add_character(i,preliminary); // add ancestral states also in descendants
		}
	    }
	    // calc branchlength
	    prefered.get_character(i) &= characters[leaf].get_character(i); // prefered is used to select along what branch a character change
	    if (prefered.get_character(i).none()) {
		++branch;
		temp=support_functions::pick_a_random_true_bit(characters[leaf].get_character(i));
		prefered.add_character(i,temp);
	    }
	}
	if (calc_branch_length) leaf->branchlength = branch;
	if (draw_ancestral_state) add_characters_as_node_comments(leaf, characters[leaf], start_char, end_char, alphabet);
	#ifdef DEBUGTREEATOR
	cerr << "Branch length:" << branch << endl;
	cerr << "Diff score: " << (characters[leaf->parent].get_score()-characters[leaf].get_score()) << endl;
	#endif //DEBUGTREEATOR
    }
    else { // For root, should not happen otherwise
	if (calc_branch_length) leaf->branchlength = 0;
	if (draw_ancestral_state) add_characters_as_node_comments(leaf, characters[leaf], start_char, end_char, alphabet);
    }
    if (leaf->left!=0) fitch_parsimony_second_pass(leaf->left, characters, prefered, calc_branch_length, draw_ancestral_state, reconstruct_tip_state, start_char, end_char, alphabet);
    if (leaf->right!=0) fitch_parsimony_second_pass(leaf->right, characters, prefered, calc_branch_length, draw_ancestral_state, reconstruct_tip_state, start_char, end_char, alphabet);
}

void tree::add_characters_as_node_comments(node* leaf, parsimony_character_vector& characters, const unsigned int start_char, const unsigned int end_char, map<char, bitset<SIZE> >& alphabet) {
    string labelcomment;
    if (leaf->nodelabel!=0) labelcomment = *leaf->nodelabel;
    labelcomment+="[&";
    //unsigned int length =characters.n_char();
    for (unsigned int i=start_char; i<=end_char; ++i) {
	if (i>0) labelcomment+=',';
	labelcomment+="trait_";
	stringstream converter;
	converter << i;
	labelcomment+=converter.str();
	labelcomment+='=';
	if (alphabet.empty())
	    labelcomment+=characters.get_character(i).to_string<char,std::string::traits_type,std::string::allocator_type>();
	else labelcomment+= alphabet::translate_bitset(characters.get_character(i),alphabet);
    }
    labelcomment+=']';
    leaf->nodelabel = nodelabels.add_string(labelcomment);
}

unsigned int tree::assign_branch_number_to_internal_nodes (node* leaf, unsigned int number) {
    if (leaf==0 || (leaf->left == 0 && leaf->right == 0)) return number;
    std::stringstream converter;
    converter << number;
    leaf->nodelabel = nodelabels.add_string(converter.str());
    ++number;
    if (leaf->left != 0) number = assign_branch_number_to_internal_nodes(leaf->left, number);
    if (leaf->right != 0) number = assign_branch_number_to_internal_nodes(leaf->right, number);
    return number;
}

bool tree::add_node_to_branch(node* leaf, node* new_node) {
    if (leaf == root || leaf == 0 || leaf->parent == 0) return false;
    node* new_internal = new node;
    new_internal->branchlength = leaf->branchlength/2;
    leaf->branchlength /= 2;
    new_internal->parent = leaf->parent;
    leaf->parent = new_internal;
    new_node->parent = new_internal;
    new_internal->left = leaf;
    new_internal->right = new_node;
    if (new_internal->parent->left == leaf) new_internal->parent->left = new_internal;
    else if (new_internal->parent->right == leaf) new_internal->parent->right = new_internal;
    else return false;
    return true;
}

void tree::clear_node_states_from_tip ( node* leaf, map<node*, parsimony_character_vector>& node_states ) {
    node* present = leaf->parent;
    while (present != 0) {
	node_states.erase(present);
	present = present->parent;
    }
}

unsigned int tree::get_tip_by_number(node* leaf, node*& tip, unsigned int tip_no) {
    if (leaf != 0) {
	if (leaf->left == 0 && leaf->right == 0) {
	    --tip_no;
	    if (tip_no == 0) tip=leaf;
	}
	else {
	    if (leaf->left!=0) tip_no = get_tip_by_number(leaf->left,tip,tip_no);
	    if (tip_no > 0 && leaf->right!=0) tip_no = get_tip_by_number(leaf->right,tip,tip_no);
	}
    }
    return tip_no;
}

void tree::stepwise_addition (vector<character_vector>& characters) {
    vector<node*> tips;
    map<node*, parsimony_character_vector > decoded_characters;
    unsigned int max_n_char = 0;
    // Map characters
    for (vector<character_vector>::iterator i= characters.begin(); i != characters.end(); ++i) {
	node* taxa = new node;
	taxa->nodelabel = nodelabels.add_string(i->get_taxon());
	tips.push_back(taxa);
	decoded_characters[taxa] = parsimony_character_vector(*i);
	unsigned int n_char = i->n_char();
	if (max_n_char != 0 && max_n_char != n_char) std::cerr << "WARNING!!! Different number of characters between taxa." << std::endl;
	if (n_char > max_n_char) max_n_char = n_char;
    }
    if (max_n_char >0) --max_n_char;
    if (!tips.empty()) {
	root->left = tips.back();
	tips.pop_back();
	root->left->parent = root;
    }
    if (!tips.empty()) {
	root->right = tips.back();
        tips.pop_back();
	root->right->parent = root;
    }
    if (!tips.empty()) {
	if (!add_node_to_branch(root->left,tips.back())) {
	    std::cerr << "Failed to add " << *tips.back()->nodelabel << " to start triplet." << std::endl;
	    delete tips.back();
	}
	tips.pop_back();
    }
    if (n_tips() < 3)
	std::cerr << "Failed to build start triplet." << std::endl;
    while (!tips.empty()) {
	if (!add_tip_to_branch_parsimony(tips.back(), decoded_characters, 0, max_n_char)) {
	    std::cerr << "Failed to add " << *tips.back()->nodelabel << " to tree." << std::endl;
	    delete tips.back();
	}
	#ifdef EXTRA
	else std::cerr << "Added " << n_tips() << " to tree." << std::endl;
	#endif // EXTRA
	//recalc_fitch_parsimony_given_added_tip(tips.back(), decoded_characters, 0, max_n_char);
	tips.pop_back();
    }
}

/*void tree::recalc_fitch_parsimony_given_added_tip ( node* leaf, map<node*, parsimony_character_vector>& node_states, unsigned int start_char, unsigned int end_char ) {
    node* prev_leaf = leaf;
    if (leaf !=0) leaf = leaf->parent;
    unsigned int score = 0;
    while (leaf != 0) {
        if (prev_leaf!= 0 && leaf->left == prev_leaf && leaf->right != 0) { 
            if (node_states.find(leaf->right) != node_states.end()) score += node_states[leaf->right].get_score();
            else score += fitch_parsimony(leaf->right, node_states, start_char, end_char);
        }
        else if (prev_leaf!= 0 && leaf->right == prev_leaf && leaf->left != 0) { 
            if (node_states.find(leaf->left) != node_states.end()) score += node_states[leaf->left].get_score();
            else score += fitch_parsimony(leaf->left, node_states, start_char, end_char);
        }
	if (node_states.find(leaf) != node_states.end()) {
	    node_states[leaf].reset_char();
	    node_states[leaf].set_score(0);
	}
	else node_states[leaf] = parsimony_character_vector();
        for (unsigned int i=start_char; i <=end_char; ++i) {
            if (leaf->left != 0 && node_states.find(leaf->left) != node_states.end() && leaf->right != 0 && node_states.find(leaf->right) != node_states.end()) { // if both descendants are initiated
                bitset<SIZE> intersect = node_states[leaf->left].get_character(i);
                intersect &= node_states[leaf->right].get_character(i);
                if (intersect.any()) node_states[leaf].add_character(i,intersect); // set to intersect
                else {
                    intersect = node_states[leaf->left].get_character(i);
                    intersect |= node_states[leaf->right].get_character(i);
                    node_states[leaf].add_character(i,intersect); // if no intersect set to union
                    ++score; // add to score
                }
            }
            else if (leaf->left != 0 && node_states.find(leaf->left) != node_states.end()) node_states[leaf].add_character(i,node_states[leaf->left].get_character(i));
            else if (leaf->right != 0 && node_states.find(leaf->right) != node_states.end()) node_states[leaf].add_character(i,node_states[leaf->right].get_character(i));
            else node_states[leaf].get_character(i).flip();
        }
	node_states[leaf].set_score(score);
        prev_leaf = leaf;
        leaf = leaf->parent;
    }
}*/

bool tree::add_tip_to_branch_parsimony(node* new_taxon, map<node*, parsimony_character_vector>& node_states, unsigned int start_char, unsigned int end_char) {
    const unsigned int number_of_tips = n_tips();
    unsigned int min_score = std::numeric_limits<unsigned int>::max();
    unsigned int best_branch = 0;
    map<node*, parsimony_character_vector> temp_nodes = node_states; // this should be unnecessary
    fitch_parsimony(root,temp_nodes,start_char,end_char);
    for (unsigned int i=1; i <= 2*number_of_tips-3; ++i) {
	unsigned int score = parsimony_score_if_tip_added_to_branch(i, temp_nodes, temp_nodes[new_taxon], start_char, end_char);
	if (score < min_score) {
	    min_score = score;
	    best_branch = i;
	}
    }
    node* branch = 0;
    if (best_branch <= number_of_tips-3) branch = get_internal_branch_by_number_unrooted(best_branch);
    else if (best_branch <= 2*number_of_tips-3) get_tip_by_number(root, branch, best_branch - (number_of_tips-3));
    else return false;
    return add_node_to_branch(branch, new_taxon);
}

unsigned int tree::parsimony_score_if_tip_added_to_branch (const unsigned int branch_no, map<node*, parsimony_character_vector>& node_states, parsimony_character_vector taxa_state, unsigned int start_char, unsigned int end_char) {
    node* branch = 0;
    const unsigned int number_of_tips = n_tips();
    if (branch_no <= number_of_tips-3) branch = get_internal_branch_by_number_unrooted(branch_no);
    else if (branch_no <= 2*number_of_tips-3) get_tip_by_number(root, branch, branch_no - (number_of_tips-3));
    else std::cerr << "Branch number larger than number of branches (" << branch_no << " > " << 2*number_of_tips-3 << ")." << std::endl;
    unsigned int score=0;
    node* prev_branch = 0;
    if (node_states.find(branch) != node_states.end()) score = node_states[branch].get_score();
    else score = fitch_parsimony(branch,node_states,start_char,end_char);
    #ifdef DEBUG
    std::cerr << "Score at start = " << score << std::endl;
    #endif //DEBUG
    while (branch != 0) {
	if (prev_branch!= 0 && branch->left == prev_branch && branch->right != 0) { 
	    if (node_states.find(branch->right) != node_states.end()) score += node_states[branch->right].get_score();
    	    else score += fitch_parsimony(branch->right, node_states, start_char, end_char);
	}
	else if (prev_branch!= 0 && branch->right == prev_branch && branch->left != 0) { 
	    if (node_states.find(branch->left) != node_states.end()) score += node_states[branch->left].get_score();
	    else score += fitch_parsimony(branch->left, node_states, start_char, end_char);
	}
	for (unsigned int i=start_char; i <=end_char; ++i) {
	    bitset<SIZE> other;
	    if (prev_branch == 0) other = node_states[branch].get_character(i);
	    else if (branch->left == prev_branch && branch->right != 0) 
		other = node_states[branch->right].get_character(i);
	    else if (branch->right == prev_branch && branch->left != 0)
		other = node_states[branch->left].get_character(i);
	    else other.flip();
	    bitset<SIZE> intersect = other;
	    intersect &= taxa_state.get_character(i);
	    if (!intersect.any()) {
		intersect = other;
		intersect |= taxa_state.get_character(i);
		++score;
	    }
	    taxa_state.add_character(i,intersect);
	}
	prev_branch = branch;
    	branch = branch->parent;
    }
    #ifdef DEBUG
    std::cerr << "Score at end = " << score << std::endl;
    #endif //DEBUG
    return score;
}

double tree::log_clade_credibility( node* leaf ) {
    double credibility(0.0);
    if (leaf != 0) {
	if (leaf->nodelabel != 0 && leaf->left != 0 && leaf->right != 0) {
	    credibility = atof(leaf->nodelabel->c_str());
	    if (credibility != 0.0) credibility = log(credibility);
	}
	if (leaf->left != 0) credibility += log_clade_credibility(leaf->left);
	if (leaf->right != 0) credibility += log_clade_credibility(leaf->right);
    }
    #ifdef DEBUG
    cerr << credibility << endl;
    #endif //DEBUG
    return credibility;
}

void tree::divide_internal_nodes_by( node* leaf, const float value ) {
    if (leaf != 0) {
	if (leaf->left != 0) divide_internal_nodes_by (leaf->left, value);
	if (leaf->right != 0) divide_internal_nodes_by (leaf->right, value);
	if (leaf->left != 0 && leaf->right != 0 && leaf->nodelabel != 0) {
	    float support = atof(leaf->nodelabel->c_str());
	    support /= value;
	    stringstream converter;
	    converter << support;
	    leaf->nodelabel = nodelabels.add_string(converter.str());
	}
    }
} 

void tree::null_short_branches( node* leaf, const double value ) {
    if (leaf != 0) {
	if (leaf->left != 0) null_short_branches( leaf->left, value);
	if (leaf->right != 0) null_short_branches( leaf->right, value);
	if (leaf->branchlength < value) leaf->branchlength = 0.0;
    }
}
// add the distance from the root of each internal node to given vector and return longest distance to a tip
double tree::get_node_heights (node* leaf, const double parent_height, vector<double>& heights) { 
    if (leaf->left != 0 || leaf->right != 0) {
	heights.push_back(parent_height+leaf->branchlength);
	double dist_left(0.0);
	double dist_right(0.0);
	if (leaf->left != 0) dist_left = get_node_heights(leaf->left, parent_height+leaf->branchlength, heights);
	if (leaf->right != 0) dist_right = get_node_heights(leaf->right, parent_height+leaf->branchlength, heights);
	if (dist_left > dist_right) return dist_left;
	else return dist_right;
    }
    else return parent_height+leaf->branchlength;
}

double tree::gamma (node* leaf) {
    vector<double> node_heights;
    node_heights.push_back(get_node_heights(leaf, 0, node_heights));
    sort(node_heights.begin(),node_heights.end());
    #ifdef DEBUG
    cerr << "N node heights: " << node_heights.size() << endl;
    #endif //DEBUG
    double sum(0.0);
    double T(0.0);
    for (unsigned int i=1; i < node_heights.size(); ++i) {
	T += (i+1) * (node_heights[i] - node_heights[i-1]);
 	if (i < node_heights.size()-1) {
	    for (unsigned int k=1; k<=i; ++k) {
		sum += (k+1) * (node_heights[k] - node_heights[k-1]);
	    }
	}
	#ifdef DEBUG
	cerr << "T: " << T << " sum: " << sum << endl;
	#endif //DEBUG
    }
    double value = (1.0/(node_heights.size()-2)) * sum;
    #ifdef DEBUG
    cerr << "value: " << value << endl;
    #endif //DEBUG
    value -= T/2.0;
    #ifdef DEBUG
    cerr << "value: " << value << endl;
    #endif //DEBUG
    value /= T * sqrt(1.0/(12*(node_heights.size()-2)));
    #ifdef DEBUG
    cerr << "value: " << value << endl;
    #endif //DEBUG
    return value;
}
/*/ functions for rate mods
bool tree::adid_rate_for_time (unsigned int rate_class, float time) {
    if (rate_class >= rates_in_time.size()) return false;
    vector<pair<unsigned int,float>>::iterator pos;
    for (pos = rate_time_cut_offs.begin(); pos != rate_time_cut_offs.end(); ++pos) {
	if (pos->second > time) break;
    }
    rates_in_time.insert(pos,make_pair(rate_class,time));
    return true;
    // add rate to rates_in_time vector in order
}
void tree::set_time_rate (unsigned int rate_class, float rate) {
    if (rate_class >= rates_in_time.size())
	while (rate_class >= rates_in_time.size()) rates.push_back(rate);
    else rates_in_time[rate_class] = rate;
}
bool tree::set_clade_rate ( unsigned int clade, unsigned int rate_class ) {
    // add checks
    if (clade >= clades.size()) return false;
    else if (rate_class >= rates.size()) return false;
    else if (clade >= clade_rate.size()) {
	while (clade >= clade_rate.size()) clade_rate.push_back(rate_class);
    } 
    else clade_rate[clade] = rate_class;
    return true;
}; 
void tree::set_rate (unsigned int rate_class, float rate) {
    if (rate_class >= rates.size())
	while (rate_class >= rates.size()) rates.push_back(rate);
    else rate[rate_class] = rate;
}
double tree::get_time_mod_branch_length (node *leaf, const float distance_to_root) {
    vector<pair<unsigned int,float>>::const_iterator pos
    const double branch_end = distance_to_root+leaf->branchlength;
    double new_br_length(0.0);
    double previous_cut_off(0.0);
    for (pos = cut_offs.begin(); pos != cut_offs.end(); ++pos) {
	#ifdef DEBUG
	cerr << "Multiplier: " << rates_in_time[pos->first] << " Cut off: " << pos->second << " New branch length: " << new_br_length << endl;
	#endif //DEBUG
	if ( distance_to_root < pos->second && branch_end > previous_cut_off) { //
	    if (distance_to_root > previous_cut_off && branch_end < pos->second) // if entirely within interval
		new_br_length += (branch_end-distance_to_root)*rates_in_time[pos->first]; // multiply entire branch
	    else if (distance_to_root < previous_cut_off && branch_end < pos->second) // if also within previous interval
		new_br_length += (branch_end-previous_cut_off)*rates_in_time[pos->first]; // multiply what is left of branch
	    else if (distance_to_root > previous_cut_off && branch_end > pos->second) // if also in next interval
		new_br_length += (pos->second-distance_to_root)*rates_in_time[pos->first]; // multiply from start of branch until end of interval
	    else if (distance_to_root < previous_cut_off && branch_end > pos->second) // if in also in both previous and next interval
		new_br_length += (pos->second-previous_cut_off)*rates_in_time[pos->first]; // multiply the distance between the intervals
	}
	#ifdef DEBUG
	cerr << "New branch length after iteration: " << new_br_length << endl;
	#endif //DEBUG
	previous_cut_off = pos->second;
    }
    if (distance_to_root < previous_cut_off && branch_end > previous_cut_off) // if sticking out beyond last interval
	new_br_length += branch_end-previous_cut_off; // add what is left
    else if (distance_to_root > previous_cut_off) new_br_length = leaf->branchlength;
    return new_br_length;
}
*/
