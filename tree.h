/********************************************************************
Copyright (C) 2014 Martin Ryberg

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
#ifndef TREEHEADER
#define TREEHEADER

#include <iostream>
#include <string.h>
#include <map>
#include <vector>
#include <set>
#include <bitset>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <sstream>
#include <ctime>
#include "string_tree.h"
#include "character_vector.h"
#include "constants.h"
#include "matrix_parser.h"
//#include "supmar/character_vector.h"
//#include "supmar/constants.h"
//#include "support_functions.cpp"
//#include "supmar/support_functions.h"

using namespace std;

/*** Class for phylogenetic trees ***/
class tree {
    protected:
	class node;

    public:
        struct poly_line_start {
            string line_start;
            int n;
            float y;
        };
        struct node_and_distance {
            node *tip;
            double distance;
        };
        tree () { //constructor 
            root = new node; //initiate the root
            srand(time(NULL));
        };
        tree ( node *leaf) { root = leaf; srand(time(NULL)); }; //initilize the class by creating the root using given node
                                             //does not set left and right too null, which may be a problem for some functions
                                             //this alternative only make sence if the given node already point to a tree
        tree (const tree& A) {
	    root = new node;
	    copy(root,A.root, 0);
            srand(time(NULL));
	}; //can not copy tree objects (yet)
        ~tree () { destroy_tree (root); }; //delete the class starting from the node
	tree operator = (const tree& A) {
	    destroy_tree(root);
	    root = new node;
	    copy(root,A.root, 0);
	    return(*this);
	};
        //parses a newik formated file and stores it in the tree (object)
	void tree_file_parser( istream& infile, map<string,string> label_translation, bool read_extra_node_annotation );
        void tree_file_parser( istream& infile ) { map<string,string> label_translation; tree_file_parser(infile,label_translation, false); };
	bool empty() {
	    if (root->left == 0 && root->right == 0) return true;
	    else return false;
	};
        void rand_topology ( const int n  );
        void print_newick( bool include_br_length, bool int_node_label ) {
            print_newick_subtree( root, include_br_length, int_node_label );
            std::cout << ";" << endl;
        };
        void print_newick( bool include_br_length ) {
            print_newick_subtree( root, include_br_length );
            std::cout << ";" << endl;
        };
        void print_newick( ) {
            print_newick_subtree( root );
            std::cout << ";" << endl;
        };
	void print_nexus_tree_intro( ostream& output, map<string,string>& translate_taxa ) { print_nexus_tree_head(output, root, translate_taxa); };
	void print_nexus_tree_intro( map<string,string>& translate_taxa ) { print_nexus_tree_intro(cout, translate_taxa); };
	void print_tree_to_nexus( ostream& output, string& name, bool include_br_length, bool int_node_label, map<string,string>& translate_taxa ) { print_tree_to_nexus_stream( output, root, name, include_br_length, int_node_label, translate_taxa); };
	void print_tree_to_nexus( ostream& output, string& name, bool include_br_length, bool int_node_label) {
	    map<string,string> translate_taxa;
	    print_tree_to_nexus_stream( output, root, name, include_br_length, int_node_label, translate_taxa);
	};
	void print_tree_to_nexus( string& name, bool include_br_length, bool int_node_label, map<string,string>& translate_taxa ) { print_tree_to_nexus( cout, name, include_br_length, int_node_label); };
	void print_tree_to_nexus( string& name, bool include_br_length, bool int_node_label) {
	    map<string,string> translate_taxa;
	    print_tree_to_nexus( cout, name, include_br_length, int_node_label, translate_taxa);
	};
	void print_nexus_tree_end( ostream& output ) { output << "End;" << endl; }
	void print_nexus_tree_end( ) { print_nexus_tree_end(cout); }
        void print_nexus( bool include_br_length, bool int_node_label ) { print_nexus_subtree( root, include_br_length, int_node_label ); };
        void print_nexus( bool include_br_length ) { print_nexus_subtree( root, include_br_length ); };
        void print_nexus( ) { print_nexus_subtree( root ); };
	void print_svg() { print_svg (true, false, 595.0, 842.0, 2.0, 1,12,"Arial",NULL); };
	void print_svg(const string* tip_color) { print_svg (true, false, 595.0, 842.0, 2.0, 1,12,"Arial",tip_color); };
        void print_svg ( bool scalebar, bool node_lable, const float width, const float height, const float offset, const unsigned int stroke_width, const int font_size, string font, const string* tip_color);
	void print_svg_no_head ( bool scalebar, bool node_lable, const float width, const float height, const float offset, const unsigned int stroke_width, const int font_size, string font, const string* tip_color, bool for_html );
	void print_svg_autoscale (bool scalebar, bool node_lable, const unsigned int font_size, const string* tip_color );
	void print_svg_autoscale_no_head (bool scalebar, bool node_lable, const unsigned int font_size, const string* tip_color, bool for_html );
        int n_nodes() { //get number of nodes
            if (root == 0) return 0; //if the root is empty there is no values to sum up
            return n_descendants(root); //calculate the sum from the root up
        };
        int n_tips () {
            if (root == 0) return 0; //if the root is empty there is no values to sum up
            return n_sub_tips(root); //calculate the number of tips desending from the root
        };
        int n_supported (float cutoff) {
            if (root == 0) return 0;
            return n_supported(root, cutoff);
        };
        void drop_tips ( const string taxa ) {
            drop_tips (root, taxa);
        };
        string tip_names() {
            return tip_names( root );
        };
        string tip_names( string separator ) {
            return tip_names( root, separator );
        };
        void midpoint_root ( );
        void outgroup_root ( const string taxa );
        void ladderize ( ) {
            ladderize ( root, 1 );
        };
        void ladderize ( bool right ) {
            ladderize ( root, right );
        };
        string not_present ( const string taxa ) {
            return not_present ( root, taxa );
        };
        bool taxon_present ( const string taxon ) {
            return taxon_present( root, taxon );
        };
        void change_tip_names ( const string tip_pairs );
        void change_tip_name( const string tip_name, const string new_tip_name ) {
            change_tip_name( root, tip_name, new_tip_name );
        };
        void multiply_br_length ( const float multiplier ) {
            multiply_br_length_subtree ( root, multiplier );
        };
	void multiply_br_length_cut_off ( const float cut_off, const float multiplier ) {
	    multiply_br_length_cut_off_subtree (root, cut_off, multiplier);
	};
	void multiply_br_length_clades ( const vector<string> &clades, const float multiplier );
        void set_br_length (const float value) {
            set_br_length_subtree ( root, value );
        };
        double sum_br_length ( ) {
            return sum_br_length( root );
        };
        double longest_branch () {
            node* node = longest_branch(root);
            return node->branchlength;
        }
	double longest_to_tip () {
	    return longest_to_tip(root);
	};
	void adjustedMPL (const string nodes_and_ages);
	void test_clock_likness() { test_clock_likness(root); }
        void print_distance_to_root_for_all ( string value_sep, string row_sep ) {
            print_distance_to_root ( root->left, 0.0, 0, value_sep, row_sep );
            print_distance_to_root ( root->right, 0.0, 0, value_sep, row_sep );
        };
        void print_distances_to_all_for_each_taxa ( const string value_divider, const string taxon_divider );
        void print_branch_lengths ( string separator, bool print_root );
        void/*tree*/ split_tree_at_longest_branch(tree* tree) { /*return*/ split_at_longest_branch ( root, tree ); };
        void/*tree*/ split_tree_at_longest_branch_unrooted(tree* tree) { /*return*/ split_at_longest_branch_unrooted ( root, tree ); };
        void split_tree_at_midpoint (tree* tree) { split_at_midpoint ( root, tree); };
        bool is_monophyletic ( string taxa ) {
            return is_monophyletic ( root, taxa );
        };
        void print_conflict_clades_reduced ( tree* tree, const float supp, bool svg );
	void internal_nodes_stat( bool sum, bool prod, bool average, bool median, bool SD, bool n_values_above, double cut_off, bool max, bool min, bool N, bool include_root ) {
	    string type;
	    internal_nodes_stat(type, sum, prod, average, median, SD, n_values_above, cut_off, max, min, N, include_root );
	}
	void internal_nodes_stat( const string& type, bool sum, bool prod, bool average, bool median, bool SD, bool n_values_above, double cut_off, bool max, bool min, bool N, bool include_root );
	unsigned int robinson_fould( tree* B );
	unsigned int shared_tips( tree* B );
	void tips_not_shared( tree* B, const string& name_tree1, const string & name_tree2 );
	void add_to_matrix_representation(map<string,vector<char> >& matrix);
	void nni (unsigned int node_no, tree& L, tree& R); // perform NNI on given branch. If branch_no == 0 random branch will be selected
	unsigned int fitch_parsimony (vector<character_vector>& characters, bool ancestral_states, bool get_branch_lengths, map<char, bitset<SIZE> >& alphabet);
	unsigned int fitch_parsimony (vector<character_vector>& characters, bool get_branch_lengths) {
	    map<char, bitset<SIZE> > alphabet;
	    return fitch_parsimony(characters, false, get_branch_lengths, alphabet);
	};
	unsigned int fitch_parsimony (vector<character_vector>& characters) {
	    map<char, bitset<SIZE> > alphabet;
	    return fitch_parsimony(characters, false, false, alphabet);
	};
	unsigned int fitch_parsimony_ancestral_states (vector<character_vector>& characters, map<char, bitset<SIZE> >& alphabet) { return fitch_parsimony(characters, true, false, alphabet); };
	void assign_branch_number_to_internal_nodes (){
	    assign_branch_number_to_internal_nodes(root, 0);
	};
	void stepwise_addition (vector<character_vector>& characters);
	void add_to_support(tree* B);
    protected:
	class node { //store the information for each node of the tree
	    public:
	    node ():branchlength(0.0), nodelabel(0), left(0), right(0), parent(0) {};
	    double branchlength; //branchlength of branch leading to node
	    string *nodelabel; //name of node
	    node *left; //left daughter node
	    node *right; //right daughter node
	    node *parent; //the parent of the node
	    //double *other;
	};
	//struct node_double_char {node* leaf; double support; char mono; };
	struct int_double2 { int n; double a; double c; };
	// Variables
	static string_tree nodelabels;
        node *root; //stores the location of the root
	//string tree_comment;
        //deletes a part of the tree from given node
        void destroy_tree (node *leaf);
	void copy (node* B, const node* A, node* parent) ;
	bool add_node_to_branch(node* leaf, node* new_node);
        //void delete_node_array( node_array* nodes );
	int print_newick_subtree( ostream& output, node *leaf, int n, bool include_br_length, bool int_node_label, map<string,string>& translate_taxa );
	int print_newick_subtree( ostream& output, node *leaf, int n, bool include_br_length, bool int_node_label) {
	    map<string,string> translate_taxa;
	    return print_newick_subtree(output,leaf,n,include_br_length,int_node_label,translate_taxa);
	};
        int print_newick_subtree( node *leaf, int n, bool include_br_length, bool int_node_label ) {
	    return print_newick_subtree (cout, leaf, n, include_br_length, int_node_label );
	};
        int print_newick_subtree( node *leaf, bool include_br_length, bool int_node_label ) {
            return print_newick_subtree( leaf, -1, include_br_length, int_node_label );
        };
        int print_newick_subtree( node *leaf, int n, bool include_br_length ) {
            return print_newick_subtree( leaf, n, include_br_length, true);
        };
        int print_newick_subtree( node *leaf, bool include_br_length ) {
            return print_newick_subtree( leaf, -1, include_br_length, true );
        };
        int print_newick_subtree( node *leaf, int n ) {
            return print_newick_subtree( leaf, n, true, true );
        };
        void print_newick_subtree( node *leaf ) {
            print_newick_subtree( leaf, -1, true, true );
        };
        void print_nexus_subtree( ostream& output, node *leaf, bool include_br_length, bool int_node_label );
        void print_nexus_subtree( node *leaf, bool include_br_length, bool int_node_label ) {
	    print_nexus_subtree( cout, leaf, include_br_length, int_node_label );
	};
        void print_nexus_subtree( node *leaf ) {
            print_nexus_subtree( leaf, true, true );
        };
        void print_nexus_subtree( node *leaf, bool include_br_length ) {
            print_nexus_subtree( leaf, include_br_length, true );
        };
        void print_nexus_tree_head( ostream& output, node *leaf, map<string,string>& translate_taxa );
	void print_tree_to_nexus_stream( ostream& output, node* leaf, string& name, bool include_br_length, bool int_node_label, map<string,string>& translate_taxa) {
	    output << "tree " << name << " = ";
	    print_newick_subtree(output, leaf, 1, include_br_length, int_node_label, translate_taxa);
	    output << ";" << endl;
	}
        void print_svg_subtree ( node* leaf, poly_line_start* branch_start, bool node_lable, int n, float x, const float x_unit, const float y_unit, const float offset, const int stroke_width, const int font_size, string font, const string* tip_color );
        int n_descendants ( node *leaf ) { //returns the number of daughter nodes for a node +1
            if (leaf == 0) return 0; //if a empty leaf return 0
            return 1 + n_descendants(leaf->left) + n_descendants(leaf->right); //return 1 plus the number of descendants to the left and right
        };
        unsigned int n_sub_tips( node *leaf ) {
            if (leaf == 0) { return 0; }
            if (leaf->left == 0 && leaf->right == 0 ) { return 1; }
            return n_sub_tips(leaf->left) + n_sub_tips(leaf->right);
        };
	unsigned int get_internal_node_by_number(node* leaf, node*& tip, unsigned int node_no); // leaft is node 1, then node number increase by one as they occure in a newick tree
	node* get_internal_branch_by_number_unrooted(unsigned int branch_no);
	unsigned int get_tip_by_number(node* leaf, node*& tip, unsigned int tip_no); // tips are numbered from the left as they occure in the tree, a pointer to the tip will be placed in tip
        bool tip_is_in_list (const string* tip, const string* list, const char separator);
        int n_supported (node* leaf, float cutoff);
        double longest_to_tip (node *leaf);
        double shortest_to_tip (node *leaf);
        void tip_furthest_from_root (node *leaf, node_and_distance* tip_distance);
        string tip_names ( node *leaf ) {
            return tip_names ( leaf, ", " );
        };
        string tip_names ( node *leaf, string separator );
	void tip_names ( node *leaf, set<string*>& tips );
	void tip_names ( node *leaf, vector<string*>& tips );
        void print_tips ( node *leaf ) {
            print_tips(leaf, "", " ");
        };
        void print_tips ( ostream& output, node *leaf ) {
            print_tips(output, leaf, "", " ");
        };
        void print_tips ( node *leaf, string leading, string trailing ) {
            print_tips( leaf, "", -1, leading, trailing );
        };
        void print_tips ( ostream& output, node *leaf, string leading, string trailing ) {
            print_tips( output, leaf, "", -1, leading, trailing );
        };
        int print_tips ( node *leaf, string leading_n, int n, string leading, string trailing ) {
	    return print_tips ( cout, leaf, leading_n, n, leading, trailing );
	};
        int print_tips ( ostream& output, node *leaf, string leading_n, int n, string leading, string trailing );
        char drop_tips (node *leaf, const string taxa);
        void re_root ( node *leaf );
        bool is_monophyletic ( node* leaf, set<string*>& taxa );
        bool is_monophyletic ( node* leaf, string& taxa );
	char is_monophyletic_unrooted (node* leaf, set<string*>& taxa, set<string*>& ignor_taxa);
	char is_monophyletic_unrooted (node* leaf, set<string*>& taxa);
	bool is_nested_in ( const node* ancestor, const node* descendent);
        node* find_midpoint_node ( node* leaf );
        void turn_nodes ( node *leaf, bool clockwise );
        inline void turn_node_clockwise ( node *leaf );
        inline void turn_node_anticlockwise ( node *leaf );
        void ladderize ( node *leaf, bool right ) {
            if (leaf->left != 0) ladderize(leaf->left, right);
            if (leaf->right !=0) ladderize(leaf->right, right);
            if (n_sub_tips(leaf->right) > n_sub_tips(leaf->left)) {
                if (right) spin_node(leaf);
            }
            else if (n_sub_tips(leaf->right) < n_sub_tips(leaf->left)) {
                if (!right) spin_node(leaf);
            }
        }
        void spin_node ( node *leaf ) {
            node *temp;
            temp = leaf->left;
            leaf->left = leaf->right;
            leaf->right = temp;
        }
	// functions for adjusted MPL //
	void adjustedMPL(map<node*,double>& given_nodeages, unsigned int n_char);
	int_double2 calculate_mean_path_root( node* leaf, map<node*,double>& given_nodeages);
	double local_adjustedMPL(node* leaf, map<node*,double>& given_nodeages, unsigned int n_char);
	int_double2 calculate_mean_path( node* leaf, const double root_age, map<node*,double>& given_nodeages);
	void assign_branches_based_on_nodeheights(node* leaf, map<node*,double>& node_depths, map<node*,double> given_nodeages); 
	int_double2 calculate_node_depth(node* leaf, const double root_age, const double path_to_root, const double rate, map<node*,double>& given_nodeages, map<node*,double>& calculated_node_depths);
	int_double2 test_clock_likness(node* leaf);
	///////////////////////////////
        string not_present ( node *leaf, const string taxa );
        bool taxon_present( node *leaf, const string taxon );
        void change_tip_name( node *leaf, const string tip_name, const string new_tip_name);
        void multiply_br_length_subtree ( node *leaf, const float multiplier );
	void multiply_br_length_cut_off_subtree (node *leaf, const float cut_off, const float multiplier );
        void set_br_length_subtree ( node *leaf, const float value );
        double sum_br_length( node *leaf );
        double two_taxa_distance ( const string* taxon1, const string* taxon2 );
        node* most_recent_common_ancestor ( const string* taxon1, const string* taxon2 );
        node* most_recent_common_ancestor ( const string taxa ); // comma separated list of taxa
	node* most_recent_common_ancestor ( set<string*> taxa );
	node* most_recent_common_ancestor ( vector<node*>& nodes );
        node* find_taxon_tip ( node *leaf, string taxon );
        node* find_taxon_tip ( node *leaf, const string* taxon );
        void print_distance_to_root ( node* leaf, double distance, int n_nodes, string value_sep, string row_sep);
        void print_branch_lengths ( node* leaf, string separator );
        node* longest_branch ( node* leaf );
        void prune_clade ( node* leaf, tree* pruned_clade );
        void split_at_longest_branch ( node* leaf, tree* tree );
        void split_at_longest_branch_unrooted ( node* leaf, tree* tree );
        void split_at_midpoint ( node* leaf, tree* tree );
	void tips_not_present_in_set (node* leaf, set<string*>& taxa, set<string*>& output); // -"-
	unsigned int tips_present_in_set (node* leaf, set<string*>& taxa);
        string supported_not_monophyletic ( set<string*>& taxa, const float supp, set<string*>& ignor_taxa ); // returns the taxa that are in conflict with the monophyly of the given taxa with given support
        void print_conflict_clades ( node* leaf, tree* tree, const float supp, set<string*>& not_in_tree, set<string*>& not_in_comp_tree, bool svg );
	char get_conflict_nodes (node* leaf, set<string*>& taxa, set<string*>& ignor_taxa, set<node*>& conflict_nodes);
	void get_internal_node_values(node* leaf, const string& type, vector<double>& values);
	void add_value_from_label_to_vector (const string& label, const string& param, vector<double>& values);
	unsigned int splits_not_in_B(node* leaf, tree* B, set<string*>& split, const unsigned int n_shared_taxa, set<string*>& not_in_A, set<string*>& not_in_B);
	void code_clades_in_matrix(node* leaf, map<string,vector<char> >& matrix, set<string>& absent_taxa);
	void add_taxa_to_matrix(node* leaf, map<string,vector<char> >& matrix, unsigned int length);
	void interchange(node* leaf, bool left);
	unsigned int fitch_parsimony (node* leaf, map<node*, parsimony_character_vector>& characters, const unsigned int start_char, const unsigned int end_char );
	void fitch_parsimony_second_pass (node* leaf, map<node*, parsimony_character_vector > characters, parsimony_character_vector prefered, bool calc_branch_length, bool draw_ancestral_state, const unsigned int start_char, const unsigned int end_char, map<char, bitset<SIZE> >& alphabet);
	unsigned int assign_branch_number_to_internal_nodes (node* leaf, unsigned int number);
	bool add_tip_to_branch_parsimony(node* new_taxon, map<node*, parsimony_character_vector>& node_states, const unsigned int start_char, const unsigned int end_char);
	unsigned int parsimony_score_if_tip_added_to_branch (const unsigned int branch_no, map<node*, parsimony_character_vector>& node_states, parsimony_character_vector taxa_state, const unsigned int start_char, const unsigned int end_char);
	void recalc_fitch_parsimony_given_added_tip ( node* leaf, map<node*, parsimony_character_vector>& node_states, unsigned int start_char, unsigned int end_char );
	void clear_node_states_from_tip ( node* leaf, map<node*, parsimony_character_vector>& node_states );
	void add_to_support(node* leaf, tree* B, set<string*>& split);
	void add_characters_as_node_comments(node* leaf, parsimony_character_vector& characters, const unsigned int start_char, const unsigned int end_char, map<char, bitset<SIZE> >& alphabet);
};
#endif //TREEHEADER
