#include "decisiveness.h"

decisiveness::decisiveness () {
    genes = 0;
    tips = 0;
    distinguished=0;
    not_distinguished=0;
    n_decisive=0;
    n_checked=0;
}
decisiveness::decisiveness ( const string* input ) {
    genes = 0;
    tips = 0;
    distinguished=0;
    not_distinguished=0;
    n_decisive=0;
    n_checked=0;
    read_genes( input );
}

void decisiveness::read_genes ( const string* input ) {
    int length = (*input).length();
    delete_gene_array (genes);
    delete_tip_genes (tips);
    genes = new gene_array;
    genes->next = 0;
    tips = new tip_genes;
    tips->next=0;
    tips->genes=0;
    int taxon = 1;
    string annotation;
    for (int i=0; i <= length; ++i) {
        if ((*input)[i]=='|' || i==length) {
            stringstream number;
            number << taxon;
            add_genes_to_tip( number.str(), &annotation );
            annotation.clear();
            ++taxon;
        }
        else annotation+=(*input)[i];
    }
}
void decisiveness::add_genes_to_tip( const string taxon, const string* input ) {
    string annotation;
    int length = (*input).length();
    for (int i=0; i <= length; ++i) {
        if ((*input)[i]==',' || i==length) {
            if (!annotation.empty()) {
                gene_array* gene_name = add_gene_to_genes( annotation );
                tip_genes* present = add_tip_genes ( taxon );
                add_gene_to_tip (gene_name, present);
                annotation.clear();
            }
            
        }
        else if ((*input)[i]==' ') continue;
        else annotation+=(*input)[i];
    }
}
decisiveness::tip_genes* decisiveness::add_tip_genes ( const string taxon ) {
    if (tips->next==0 && tips->lable.empty()) {
        tips->lable=taxon;
        return tips;
    }
    tip_genes* present=get_tip_pointer( taxon );
    if (present==0) {
        present=tips;
        while (present->next!=0) present=present->next;
        present->next= new tip_genes;
        present=present->next;
        present->next=0;
        present->genes=0;
        present->lable=taxon;
    }
    return present;
}
void decisiveness::add_gene_to_tip (gene_array* gene, tip_genes* tip) {
    if (tip->genes==0) {
        tip->genes = new node_genes;
        tip->genes->next = 0;
        tip->genes->gene=0;
    }
    add_gene_to_node_genes (gene,tip->genes);
}
void decisiveness::add_gene_to_node_genes ( gene_array* gene, node_genes* node ) {
    if (gene != 0) {
        while (node->next !=0 && node->gene!=gene) node=node->next;
        if (node->gene==0) node->gene = gene;
        else if (node->gene!=gene) {
            node->next = new node_genes;
            node->next->next = 0;
            node->next->gene = gene;
        }
    }
}
decisiveness::tip_genes* decisiveness::get_tip_pointer ( const string taxon ) {
    tip_genes* present=tips;
    while (present->lable.compare(taxon) && present->next!=0) present=present->next;
    if (!present->lable.compare(taxon)) return present;
    else return 0;
}
decisiveness::gene_array* decisiveness::add_gene_to_genes (const string gene) {
    if (genes->name.empty() && genes->next==0) {
        genes->name=gene;
        return genes;
    }
    else {
        gene_array* present=genes;
        while (present->name.compare(gene) && present->next!=0) present=present->next;
        if (present->name.compare(gene) && present->next==0) {
            present->next=new gene_array;
            present = present->next;
            present->next=0;
            present->name = gene;
        }
        return present;
    }
}
void decisiveness::delete_gene_array ( gene_array* node) {
    if (node !=0) {
        if (node->next!=0) delete_gene_array(node->next);
        delete node;
    }
}
void decisiveness::delete_tip_genes ( tip_genes* node) {
    if (node !=0) {
        if (node->next!=0) delete_tip_genes(node->next);
        delete_node_genes (node->genes);
        delete node;
    }
}
void decisiveness::delete_node_genes ( node_genes* node ) {
    if (node !=0) {
        if (node->next != 0) delete_node_genes (node->next);
        delete node;
    }
}
decisiveness::node_genes* decisiveness::get_genes_sub_tree ( const node* leaf ) {
    if (leaf == 0) return 0;
    else if ( leaf->left == 0 && leaf->right == 0 ) {
        tip_genes* tip_genes = get_tip_pointer(*leaf->nodelabel);
        if (tip_genes != 0) return join_node_genes(tip_genes->genes,0);
        else return 0;
    }
    else {
        node_genes* left_genes = 0;
        if ( leaf->left != 0 ) left_genes = get_genes_sub_tree( leaf->left );
        node_genes* right_genes = 0;
        if ( leaf->right != 0 ) right_genes = get_genes_sub_tree( leaf->right );
        node_genes* joint = join_node_genes( left_genes, right_genes );
        delete_node_genes(left_genes);
        delete_node_genes(right_genes);
        return joint;
    }
}
decisiveness::node_genes* decisiveness::join_node_genes ( node_genes* one, node_genes* two ) {
    node_genes* joint=0;
    while (one != 0) {
        if (joint==0) {
            joint = new node_genes;
            joint->next=0;
            joint->gene=0;
        }
        add_gene_to_node_genes(one->gene,joint);
        one=one->next;
    }
    while (two != 0 ) {
        if (joint==0) {
            joint = new node_genes;
            joint->next=0;
            joint->gene=0;
        }   
        add_gene_to_node_genes(two->gene,joint);
        two = two->next;
    }
    return joint;
}
bool decisiveness::present_in_node_genes ( node_genes* genes, gene_array* gene ) {
    if (gene != 0) {
        while (genes != 0) {
            if (genes->gene == gene) return true;
            genes=genes->next;
        }
    }
    return false;
}
bool decisiveness::is_distinguished ( node_genes* one, node_genes* two, node_genes* three, node_genes* four ) {
    node_genes* intersect_one=0;
    node_genes* intersect_two=0;
    while (two!=0) {
        if (present_in_node_genes(one, two->gene)) {
            if (intersect_one==0) {
                intersect_one = new node_genes;
                intersect_one->next=0;
                intersect_one->gene=0;
            }
            add_gene_to_node_genes(two->gene,intersect_one);
        }
        two=two->next;
    }
    if (intersect_one != 0) {
        while (three!=0) {
            if (present_in_node_genes(intersect_one, three->gene)) {
                if (intersect_two==0) {
                    intersect_two = new node_genes;
                    intersect_two->next=0;
                    intersect_two->gene=0;
                }
                add_gene_to_node_genes(three->gene,intersect_two);
            }
            three=three->next;
        }
    }
    delete_node_genes(intersect_one);
    if (intersect_two != 0) {
        while (four!=0) {
            if (present_in_node_genes(intersect_two, four->gene)) {
                delete_node_genes(intersect_two);
                return true;
            }
            four=four->next;
        }
    }
    delete_node_genes(intersect_two);
    return false;
}
decisiveness::node_genes* decisiveness::sub_check ( const node* leaf, node_genes* parent_parent, node_genes* parent_other ) {
    if (leaf == 0) return 0;
    if ( leaf->left == 0 && leaf->right == 0 ) {
        tip_genes* tip_genes = get_tip_pointer(*leaf->nodelabel);
        if (tip_genes != 0) return join_node_genes(tip_genes->genes,0);
        else return 0;
    }
    else {
        node_genes* left_genes = 0;
        node_genes* right_genes = 0;
        node_genes* parents = join_node_genes(parent_parent,parent_other);
        if ( leaf->right != 0) right_genes=get_genes_sub_tree( leaf->right );
        if ( leaf->left != 0) left_genes = sub_check( leaf->left, parents, right_genes );
        if ( parent_parent != 0 &&  parent_other!= 0  && left_genes != 0 && right_genes != 0) {
            if (is_distinguished(parent_parent,parent_other,left_genes,right_genes)) ++distinguished;
            else ++not_distinguished;
        }
        
        if ( leaf->right != 0) sub_check( leaf->right, parents, left_genes );
        delete_node_genes(parents);
        node_genes* return_genes=join_node_genes(left_genes,right_genes);
        delete_node_genes(left_genes);
        delete_node_genes(right_genes);
        return return_genes;
    }
}
bool decisiveness::check_tree ( ) {
    int prev_not_distinguished=not_distinguished;
    node_genes* right = 0;
    node_genes* left=0;
    if ( root->right->right != 0) right=get_genes_sub_tree( root->right->right );
    if ( root->right->left != 0) left=get_genes_sub_tree( root->right->left );
    if (root->right->right == 0 && root->right->left == 0 ) right=get_genes_sub_tree( root->right );
    node_genes* left_nodes=0;
    if ( root->left != 0) left_nodes=sub_check( root->left, right, left );
    delete_node_genes(right);
    delete_node_genes(left);
    right=0;
    left=0;
    //if ( root->left->right != 0) right=get_genes_sub_tree( root->left->right );
    //if ( root->left->left != 0) left=get_genes_sub_tree( root->left->left );
    //if (root->left->right == 0 && root->left->left == 0 ) right=get_genes_sub_tree( root->right );
    //right=get_genes_sub_tree( root->right );
    if ( root->right != 0 ) right=sub_check( root->right, left_nodes, 0 );
    delete_node_genes(left_nodes);
    delete_node_genes(right);
    if (prev_not_distinguished-not_distinguished == 0) return true;
    else return false;
}

void decisiveness::on_random_tree( int n_trees ) {
    int n_taxa=0;
    tip_genes* present=tips;
    while (present != 0) {
        ++n_taxa;
        present = present->next;
    }
    for (int i=0; i < n_trees; ++i) {
        rand_topology(n_taxa);
        if (check_tree ( )) ++n_decisive;
        ++n_checked;
    }
}
float decisiveness::get_decisiveness ( ) {
    return float(n_decisive)/n_checked;
}
float decisiveness::get_distinguished ( ) {
    return float(distinguished)/(distinguished+not_distinguished);
}
