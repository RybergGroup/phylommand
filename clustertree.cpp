#include "clustertree.h"

bool clustertree::short_br_cluster (node *leaf, const double cut_off) {
    if (leaf->left ==0 && leaf->right ==0) return 1;
    bool left_cluster=0;
    bool right_cluster=0;
    if (leaf->left !=0) left_cluster = short_br_cluster(leaf->left, cut_off);
    if (leaf->right !=0) right_cluster = short_br_cluster(leaf->right, cut_off);
    if (left_cluster && right_cluster && (shortest_to_tip(leaf->left) + shortest_to_tip(leaf->right)) < cut_off ) return 1;
    else {
        if (left_cluster) { print_tips(leaf->left); std::cout << endl; }
        if (right_cluster) { print_tips(leaf->right); std::cout << endl; }
    }
    return 0;
}

void clustertree::br_length_clust_max_cout ( const float cut_off ) { // function that 'split tree' on branches longer than cut-off (sum off branches from root) and print tip labels for each subtree
    cluster_struct clusters;
    float rooted_br_length = 0.0;
    if ( root->left == 0 && root->right == 0) clusters.add_tip( root );
    else {
        if ( root->left != 0 ) rooted_br_length += root->left->branchlength;
        if ( root->right != 0 ) rooted_br_length += root->right->branchlength;
        if ( rooted_br_length > cut_off ) {
            if ( root->left != 0 ) {
                clusters.add_node();
                br_length_clust_max_cout ( root->left, &clusters, cut_off );
                clusters.move_to_parent();
            }
            if ( root->right != 0 ) {
                clusters.add_node();
                br_length_clust_max_cout ( root->right, &clusters, cut_off );
                clusters.move_to_parent();
            }
        }
        else {
            if ( root->left != 0 ) {
                br_length_clust_max_cout ( root->left, &clusters, cut_off );
            }
            if ( root->right != 0 ) {
                br_length_clust_max_cout ( root->right, &clusters, cut_off );
            }
        }
    }
    clusters.print_clusters();
}

void clustertree::br_length_clust_max_cout ( node* leaf, cluster_struct* clusters, const float cut_off ) { //supports above function
    if (leaf->branchlength > cut_off) clusters->add_node();
    if (leaf->left == 0 && leaf->right == 0) clusters->add_tip( leaf );
    else {
        if ( leaf->left !=0 ) br_length_clust_max_cout ( leaf->left, clusters, cut_off );
        if ( leaf->right !=0 ) br_length_clust_max_cout ( leaf->right, clusters, cut_off );
    }
    if (leaf->branchlength > cut_off) clusters->move_to_parent();
}

void clustertree::br_length_clust_max_cout ( float cut_off, int max_size ) { // fully functionel old function that is not in use anymore, can give strange results on 'unrooted' trees
    if ( (root->left->branchlength + root->right->branchlength) > cut_off ) {
        if ( n_sub_tips(root->left) <= n_sub_tips(root->right) && n_sub_tips(root->left) < max_size) {
            print_tips( root->left );
            cout << endl;
        }
        else {
            if ( root->left != 0 ) br_length_clust_max_cout( root->left, cut_off, max_size ); 
        }
        if ( n_sub_tips(root->left) >= n_sub_tips(root->right) && n_sub_tips (root->right) <= max_size ) {
            print_tips( root->right );
            cout << endl;
        }
        else {
            if ( root->right != 0 ) br_length_clust_max_cout( root->right, cut_off, max_size );
        }
    }
    else {
        if ( root->left != 0 ) br_length_clust_max_cout( root->left, cut_off, max_size );
        if ( root->right != 0 ) br_length_clust_max_cout( root->right, cut_off, max_size );
    }
}

void clustertree::br_length_clust_max_cout( node *leaf, float cut_off, int max_size ) { // support above function
    if ( leaf != 0 ) {
        if ( leaf->parent != root && leaf->branchlength > cut_off && n_sub_tips( leaf ) <= max_size ) {
            print_tips ( leaf );
            cout << endl;
        }
        else {
            if ( leaf->left != 0 ) br_length_clust_max_cout( leaf->left, cut_off, max_size );
            if ( leaf->right != 0 ) br_length_clust_max_cout( leaf->right, cut_off, max_size );
        }
    }
}

string clustertree::name_clust_cout ( unsigned int name_pos, char separator, node *leaf ) {
    if ( leaf->left == 0 && leaf->right == 0 ) {
        string species;
        unsigned int i = 1;
        unsigned int position = 0;
	unsigned int string_length = leaf->nodelabel->length();
        while ( i < name_pos && position < string_length) {
            if ( leaf->nodelabel->at(position) == separator ) i++;
            ++position;
        }
        while ( position < string_length && leaf->nodelabel->at(position) != separator ) { //leaf->nodelabel->at(position) != '\0') {
            species += leaf->nodelabel->at(position);
            ++position;
        }
        return species;
    }
    else {
        string species_left;
        if ( leaf->left != 0 ) species_left = name_clust_cout ( name_pos, separator, leaf->left );
        string species_right;
        if ( leaf->right !=0 ) species_right = name_clust_cout ( name_pos, separator, leaf->right );
        if (species_left.compare("Clade with different names") != 0 && species_right.compare("Clade with different names") != 0 && species_left.compare(species_right) == 0 ) return species_left;
        else {
            if ( leaf->left != 0 && species_left.compare("Clade with different names") != 0 ) {
                print_tips ( leaf->left );
                cout << endl;
            }
            if ( leaf->right != 0 && species_right.compare("Clade with different names") != 0 ) {
                print_tips ( leaf->right );
                cout << endl;
            }
            return "Clade with different names";
        }
    }
}

#ifdef DATABASE
string clustertree::db_clust_cout ( const char* database, string table, string key, string column, node *leaf ) {
    if ( leaf->left == 0 && leaf->right == 0 ) {
        string clust_string;
        sqlite3 *db;
        if ( sqlite3_open(database, &db) == 0 ) {
            sqlite3_stmt *statement;
            string query = "SELECT ";
            query += column;
            query += " FROM ";
            query += table;
            query += " WHERE ";
            query += key;
            query += "='";
            query += *(leaf->nodelabel);
            query += "';";
            char *c_query=new char[query.size()+1];
            c_query[query.size()]=0;
            memcpy(c_query,query.c_str(),query.size());
            
            if(sqlite3_prepare_v2(db, c_query, -1, &statement, 0) == SQLITE_OK) {
                delete[] c_query;
                if (sqlite3_step(statement) == SQLITE_ROW) {
                    const unsigned char* char_string;
                    char_string = sqlite3_column_text(statement,0);
                    int i=0;
                    while (char_string[i] != 0) {
                        clust_string += char_string[i];
                        ++i;
                    }
                }
                else {
                    std::cerr << "WARNING!!! Could not retrieve database entry (out of bound)." << endl;
                    return "Clade with different names";
                }
                sqlite3_finalize(statement);
            }
            else {
                std::cerr << "WARNING!!! Could not retrieve database entry (query error)." << endl;
                delete[] c_query;
                return "Clade with different names";
            }
        }
        else {
            std::cerr << "WARNING!!! Could not open database " << database << "." << endl;
            return "Clade with different names";
        }
        sqlite3_close(db);
        return clust_string;
    }
    else {
        string clust_string_left;
        if ( leaf->left != 0 ) clust_string_left = db_clust_cout ( database, table, key, column, leaf->left );
        string clust_string_right;
        if ( leaf->right !=0 ) clust_string_right = db_clust_cout ( database, table, key, column, leaf->right );
        if (clust_string_left.compare("Clade with different names") != 0 && clust_string_right.compare("Clade with different names") != 0 && clust_string_left.compare(clust_string_right) == 0 ) return clust_string_left;
        else {
            if ( leaf->left != 0 && clust_string_left.compare("Clade with different names") != 0 ) {
                print_tips ( leaf->left );
                std::cout << "[" << clust_string_left << "]" << endl;
            }
            if ( leaf->right != 0 && clust_string_right.compare("Clade with different names") != 0 ) {
                print_tips ( leaf->right );
                std::cout << "[" << clust_string_right << "]" << endl;
            }
            return "Clade with different names";
        }
    }
}
#endif /*DATABASE*/

/*********************************************/
/*** Functions of the cluster_struct class ***/
/*********************************************/
void clustertree::cluster_struct::delete_struct (cluster* node) { // function to delete all allocated memory, used by destructor
    delete_list (node->tipnames);
    delete_children (node->branches);
    delete node;
}

clustertree::cluster_struct::cluster_struct () { // constructor
    root = new cluster;
    present = root;
    root->parent = 0;
    root->tipnames=0;
    root->branches = 0;
}

void clustertree::cluster_struct::delete_list ( list* node ) {
    if (node != 0) {
        if (node->next != 0) delete_list (node->next);
        delete node;
    }
}
void clustertree::cluster_struct::delete_children ( children* node ) {
    if (node != 0) {
        if (node->next) delete_children (node->next);
        if (node->child !=0) delete_struct(node->child);
        delete node;
    }
}

void clustertree::cluster_struct::add_tip ( tree::node* name) { // add a node to the tipnames list of present cluster
    if (present->tipnames == 0) {
        present->tipnames = new list;
        present->tipnames->name = name;
        present->tipnames->next = 0;
    }
    else {
        list*  enter = present->tipnames;
        while (enter->next != 0) enter = enter->next;
        enter->next = new list;
        enter->next->name = name;
        enter->next->next = 0;
    }
}

void clustertree::cluster_struct::add_node ( ) { // create a new node and set present to be that node
    if (present->branches == 0) {
        present->branches = new children;
        present->branches->child = new cluster;
        present->branches->child->parent = present;
        present->branches->child->tipnames = 0;
        present->branches->child->branches = 0;
        present->branches->next = 0;
        present = present->branches->child;
    }
    else {
        children* enter = present->branches;
        while (enter->next != 0) enter = enter->next;
        enter->next = new children;
        enter->next->child = new cluster;
        enter->next->child->parent = present;
        enter->next->child->tipnames = 0;
        enter->next->child->branches = 0;
        enter->next->next = 0;
        present = enter->next->child;
    }
}
void clustertree::cluster_struct::print_clusters ( cluster* node ) { // print clusters, each tip separated by space each cluster by newline
    if (node->tipnames != 0) {
        list* enter = node->tipnames;
        while (enter->next != 0) {
            if (enter->name != 0) std::cout << *(enter->name->nodelabel) << " ";
            enter=enter->next;
        }
        if (enter->name != 0) std::cout << *(enter->name->nodelabel) << endl;
    }
    if (node->branches != 0) {
        children* enter = node->branches;
        while (enter->next != 0) {
            if (enter->child != 0) print_clusters (enter->child);
            enter=enter->next;
        }
        if (enter->child != 0) print_clusters (enter->child);
    }
}

