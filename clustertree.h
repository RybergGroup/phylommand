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

#include "tree.h"
#ifdef DATABASE
#include "sqlite3.h"
#endif /*DATABASE*/
using namespace std;

class clustertree : public tree {
    public:
    class cluster_struct { // class to store tree nodes in clusters arranged in a pollyfurcating tree
	public:
	    //struct cluster;
	    struct cluster; // Forward declare to be able to declare children
	    struct list { // to store several tree nodes per cluster
		tree::node* name;
		list* next;
	    };
	    struct children { // to store clusters some how related to the present cluster
		cluster* child;
		children* next;
	    };
	    struct cluster {
		cluster* parent;
		list* tipnames;
		children* branches;
	    };
	    cluster_struct ();
	    ~cluster_struct ( ) { delete_struct(root); };
	    void add_tip ( tree::node* name); // add a node to the tipnames list of present cluster
	    void add_node ( ); // create a new node and set present to be that node
	    void move_to_parent () { present = present->parent; };
	    void print_clusters ( ) { print_clusters ( root ); }
	private:
	    cluster* root;
	    cluster* present;
	    void delete_struct (cluster* node);
	    void delete_list ( list* node );
    	    void delete_children ( children* node );
	    void print_clusters ( cluster* node );
    };

        void short_br_clust ( const double cut_off ) {
            if (short_br_cluster (root, cut_off)) { print_tips(root); std::cout << endl; }
        };
        void br_length_clust_max_cout( const float cut_off, const unsigned int max_size );
        void br_length_clust_max_cout( float cut_off );
        void name_clust_cout ( unsigned int name_pos, char separator ) {
            name_clust_cout( name_pos, separator, root );
        };
        void name_clust_cout ( ) {
            name_clust_cout( 0, ' ', root );
        };
        #ifdef DATABASE
        void db_clust_cout ( const char* database, string table, string key, string column ) {
            string cluster = db_clust_cout ( database, table, key, column, root );
            if (cluster.compare("Clade with different names")!=0) {
                if ( root->left != 0 || root->right !=0 ) {
                    print_tips ( root->left );
                    print_tips ( root->right );
                }
                else std::cout << *root->nodelabel << ' ';
                std::cout << "[" << cluster << "]" << endl;
            }
        }
        #endif /*DATABASE*/
    private:
        void br_length_clust_max_cout( node *leaf, const float cut_off, const unsigned int max_size );
        void br_length_clust_max_cout ( node* leaf, cluster_struct* clusters, const float cut_off );
        string name_clust_cout ( unsigned int name_pos, char separator, node *leaf ); 
        #ifdef DATABASE
        string db_clust_cout ( const char* database, string table, string key, string column, node *leaf );
        #endif /*DATABASE*/
        bool short_br_cluster (node *leaf, const double cut_off);
};
