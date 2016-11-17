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
#ifndef NJTREEHEADER
#define NJTREEHEADER

#include "tree.h"
using namespace std;

class njtree : public tree {
    public:
        /*njtree( ) {
            //distance_matrix=0;
        };
        ~njtree( ) {
            delete_node_and_distance_array( distance_matrix );
        };*/
        struct node_and_distance {
            vector<float> distances;
            float S;
            node *child;
        };
	unsigned int n_nodes_in_array() { return distance_matrix.size(); }
        void read_distance_matrix( istream& infile, bool lables);
        void build_nj_tree (bool quiet);
        void build_nj_tree () { build_nj_tree(true); };
        void print_node_and_distance_array() { print_node_and_distance_array ( distance_matrix.begin(), cout ); };
        void print_node_and_distance_array( ostream& output) {
            print_node_and_distance_array ( distance_matrix.begin(), output );
        }
	bool matrix_good();
    private:
        vector<node_and_distance> distance_matrix;
        void print_node_and_distance_array ( vector<node_and_distance>::const_iterator start_node, ostream& output );
	void print_unfinished_tree( ostream& output );
};

#endif //NJTREEHEADER
