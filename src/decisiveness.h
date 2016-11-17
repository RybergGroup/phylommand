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
using namespace std;

class decisiveness : public tree {
    public:
    struct gene_array {
        string name;
        gene_array* next;
    };
    struct node_genes {
        gene_array* gene;
        node_genes* next;
    };
    struct tip_genes {
        string lable;
        node_genes* genes;
        tip_genes* next;
    };
    decisiveness ();
    decisiveness ( const string* input );
    ~decisiveness () {
        delete_gene_array (genes);
        delete_tip_genes (tips);
    };
    void read_genes ( const string* input );
    void on_random_tree( int n_trees );
    float get_decisiveness ( );
    float get_distinguished ( );
    private:
    int distinguished; // the number of nodes that are been deistinguished
    int not_distinguished; // the number of nodes that are not distinguished
    int n_decisive; // number of trees that are decissive
    int n_checked;
    gene_array* genes; // the names of the genes
    tip_genes* tips; // the genes for each tip
    void add_genes_to_tip( const string taxon, const string* input ); // add a comma separated string of genes to 
    tip_genes* add_tip_genes ( const string taxon ); // add a tip to the array of tips
    void add_gene_to_tip ( gene_array* gene, tip_genes* tip); // add a gene to the specified tip
    void add_gene_to_node_genes ( gene_array* gene, node_genes* node ); // add gene to a node_genes array
    tip_genes* get_tip_pointer ( const string taxon ); // get the pointer to the tip with given name
    gene_array* add_gene_to_genes (const string gene); // add a gene to genes array
    void delete_gene_array ( gene_array* node); // delete the gene array
    void delete_tip_genes ( tip_genes* node); // delete the tips array
    void delete_node_genes ( node_genes* node ); // delete the gene array of a node
    node_genes* get_genes_sub_tree ( const node* leaf );
    node_genes* join_node_genes ( node_genes* one, node_genes* two );
    bool present_in_node_genes ( node_genes* genes, gene_array* gene );
    bool is_distinguished ( node_genes* one, node_genes* two, node_genes* three, node_genes* four );
    node_genes* sub_check ( const node* leaf, node_genes* parent_parent, node_genes* parent_other );
    bool check_tree ( );
};
