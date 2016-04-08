/********************************************************************
Copyright (C) 2011 Martin Ryberg

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

contact: kryberg@utk.edu
*********************************************************************/
#ifndef NJTREEHEADER
#define NJTREEHEADER

#include "tree.h"
using namespace std;

class njtree : public tree {
    public:
        njtree( ) {
            distance_matrix=0;
        };
        ~njtree( ) {
            delete_node_and_distance_array( distance_matrix );
        };
        struct distance_array {
            distance_array *next;
            float distance;
        };
        struct node_and_distance_array {
            distance_array *distances;
            float S;
            node *child;
            node_and_distance_array *next;
        };
        void read_distance_matrix( istream& infile, bool lables);
        void build_nj_tree ();
        void print_node_and_distance_array() {
            print_node_and_distance_array ( distance_matrix );
        }
    private:
        node_and_distance_array* distance_matrix;
        int n_taxa_node_and_distance_array ( node_and_distance_array* array);
        void delete_node_and_distance_array(node_and_distance_array* start_node);
        void delete_node_and_distance_array_node(node_and_distance_array* delete_node, node_and_distance_array* previous_node);
        void delete_distance_array_node ( distance_array* delete_node, distance_array* previous_node, node_and_distance_array* father );
        void delete_distance_array ( distance_array* start_node );
        void print_node_and_distance_array ( node_and_distance_array* start_node );
        void print_distance_array ( distance_array* start_node );
        void add_to_distance_array ( float value, distance_array* start_node );
};

#endif //NJTREEHEADER
