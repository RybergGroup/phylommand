#include <iostream>
#include <stdlib.h>
#include <string.h>
//#include "sqlite3.h"

using namespace std;

template<typename kind>
struct list {
    kind *daughter;
    list *next;
};

const int precision = 1000; // The precision of the value stored: 1/precision
const int heighest = 2;     // The highest value that can be stored, higher values will be rounded to this value
const int array_length = precision * heighest;
const float cutoff = 0.01;  // The MAD cut off for what is considered alignable
//const int max_group_size = 1000;
const char gb_data[] = {'g','b','_','d','a','t','a','\0'}; // The name of the table with the original GenBank data

class align_group {
    private:
        struct node {
            list<node> *nodes;
            node *parent;
            string taxon;
            int saturation_values[array_length];
        };
        node *root;
    public:
        align_group ( ){
            root = new node;
            root->nodes = 0;
            for (int i=0; i < array_length; ++i) root->saturation_values[i] = 0;
        };
        ~align_group ( ) {
            destroy_sub_tree( root );
        };
        void insert_value ( const string taxon_string1, const string taxon_string2, float value );
/*        void insert_value ( string accno1, string accno2, float value, sqlite3 *db ) {
            string taxon_string1=get_taxon_string(accno1, db);
            string taxon_string2=get_taxon_string(accno2, db);
            insert_value ( taxon_string1, taxon_string2, value );
        };*/

        string get_levels ( string gene ) {
            return get_levels ( root, gene); //db, gene );
        };
        float aprox_mad ( ) {
            int values[array_length];
            for (int i=0; i < array_length; ++i) values[i]=0;
            add_values (values, root);
            return calc_aprox_mad( values );
        };
        void print_hierarchy ( ) {
            print_hierarchy ( root );
        };
    private:
        void destroy_sub_tree ( node *leaf ) {
            if (leaf->nodes != 0) destroy_sub_list ( leaf->nodes );
            delete leaf;
        };
        void destroy_sub_list ( list<node> *entry ) {
            if ( entry != 0 ) {
                destroy_sub_list( entry->next );
                destroy_sub_tree (entry->daughter);
                delete entry;
            }
        };
        void find_node_insert_value ( string taxon, float value, node *leaf );
        void add_node ( list<node> *position, node *parent, string taxon );
        string get_levels ( node *leaf, string gene); //sqlite3 *db, string gene );
        void add_values ( int values[], node *leaf );
        float calc_aprox_mad( int values[] );
        //int get_n_taxa( string taxon, sqlite3 *db, string table );
        // function for database handling
        //string get_taxon_string( string accno, sqlite3 *db );
        // Function to print the values for each level in the hierarchy, mostly for debugging purposes
        void print_hierarchy ( node *leaf );
};
