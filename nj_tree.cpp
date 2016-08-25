#include "nj_tree.h"

int njtree::n_taxa_node_and_distance_array ( node_and_distance_array* array) {
    int n=1;
    while (array->next!=0) {
        ++n;
        array=array->next;
    }
    return n;
}

void njtree::build_nj_tree ( ) {
    #ifdef DEBUG
    print_node_and_distance_array(cerr);
    #endif //DEBUG
    int r = n_taxa_node_and_distance_array(distance_matrix);
    while (r>2) {
        node_and_distance_array* present=distance_matrix;
        int n=0;
	#ifdef DEBUG
        cerr << "Step 1: calculate sum of distances for each taxa" << endl;
	#endif //DEBUG
	//#ifdef DEBUG
	//print_node_and_distance_array(cerr);
	//#endif //DEBUG
        while (present!=0) {
	    #ifdef DEBUG
	    if (present != 0 && present->child != 0)
		cerr << "calculating for: " << *(present->child->nodelabel) << " (" << present << ";" << present->next << ")" << endl;
	    else cerr << "No matrix." << endl;
	    #endif //DEBUG
	    node_and_distance_array* previous = distance_matrix;
	    present->S = 0;
	    int n_taxa=1;
	    // Before reaching row for the taxa the distances
	    // are summarized for, find the distance to the taxa
	    // for which row we are in
	    while (previous != present ) {
		distance_array* cell=previous->distances;
		for (int i=0; i<n-n_taxa; ++i) cell=cell->next;
		present->S += cell->distance;
		previous=previous->next;
		++n_taxa;
	    }
	    // Once we reach the row for the taxa the distances
	    // are summarized for the rest of the distances are
	    //  in that row
	    distance_array* cell=present->distances;
	    while (cell!=0 && n_taxa-n > 0) {
		present->S += cell->distance;
		cell=cell->next;
	    }
	    ++n;
	    present=present->next;
	}
	#ifdef DEBUG
        cerr <<  "Step 2: calculate pair with smallest M" << endl;
	#endif //DEBUG
        float M=0; // values in Q matrix
        int i=0;
        int j=0;
        n=0;
        present=distance_matrix;
	//#ifdef DEBUG
	//print_node_and_distance_array(cerr);
	//#endif //DEBUG
        while (present != 0) {
	    #ifdef DEBUG
	    cerr << "calculating for: " << *(present->child->nodelabel) << " (" << present << ";" << present->next << ")" << endl;
	    print_distance_array( present->distances, cerr );
	    cerr << endl;
	    #endif //DEBUG
            distance_array* cell=present->distances;
            int n_taxa=0;
            while (cell!=0 && n_taxa-n > 0) {
                //float Sj;
                node_and_distance_array* j_node = present->next;
                for (int x=0; x<n_taxa; ++x) j_node = j_node->next;
                float value=(r-2)*cell->distance - present->S - j_node->S;
		#ifdef DEBUG
		cerr << value << " (" << cell->distance << ";" << cell << ";" << cell->next << ") ";
		#endif //DEBUG
                if (value < M) { // No need to save all values just find lowest
                    M=value;
                    i=n;         // note the first
                    j=n_taxa;    // and second taxa in the pair
                }
                ++n_taxa;
                cell=cell->next;
            }
    	    #ifdef DEBUG
	    cerr << endl;
	    #endif //DEBUG
            ++n;
            present=present->next;
        }
	#ifdef DEBUG
	cerr << endl;
        cerr << "Step 3-4: Create a node joining selected taxa" << endl;
	#endif //DEBUG
        node* new_node = new node;
        node_and_distance_array* join_node1 = distance_matrix;
        for (int x=0; x<i; ++x) join_node1 = join_node1->next;
        new_node->left = join_node1->child;
        join_node1->child=0;
        node_and_distance_array* join_node2 = join_node1->next;
        for (int x=0; x<j; ++x) join_node2 = join_node2->next;
	#ifdef DEBUG
	if (join_node1->child != 0 && join_node2->child != 0) 
	    cerr << "Creating node to join " <<  join_node1->child->nodelabel << " and " << join_node2->child->nodelabel << endl;
	#endif //DEBUG
        new_node->right = join_node2->child;
        join_node2->child=0;
        float length;
        distance_array* cell=join_node1->distances;
        for (int x=0; x<j; ++x) cell = cell->next;
        length=cell->distance;
        new_node->left->branchlength= (length/2) + (join_node1->S - join_node2->S)/(2*(r-2));
        new_node->right->branchlength=length-new_node->left->branchlength;
        new_node->left->parent=new_node;
        new_node->right->parent=new_node;
	#ifdef DEBUG
        cerr << "Step 5: calculate new distance matrix" << endl;
	#endif //DEBUG
        // create a new row
        node_and_distance_array* new_distance_node = new node_and_distance_array;
        new_distance_node->child=new_node; // point it to new node
        new_distance_node->S=0;
        new_distance_node->next=distance_matrix; // set it before all other rows
        new_distance_node->distances = new distance_array; // give it a distance array
        new_distance_node->distances->next=0;
        new_distance_node->distances->distance=0;
        distance_matrix = new_distance_node;
        /**/
        distance_array* distant_taxa=0;
        distance_array* cell_i = 0;
        distance_array* cell_j = 0;
        present = distance_matrix->next;
        for (int taxa=0; taxa<r; ++taxa) { // for each taxa
            if (taxa < i) { // if taxa before i the distance between them is in taxas row
                cell_i=present->distances;
                for (int x=0; x <i-taxa-1; ++x)  cell_i=cell_i->next;
            }
            else { // otherwise it is in its row
                if (taxa==i) {
                    cell_i=present->distances;
                    present=present->next;
                    continue;
                }
                else if (taxa==i+1);
                else cell_i=cell_i->next; // step one cell furter in the row
            }
            if (taxa < i+j+1) { // same for j
                cell_j=present->distances;
                for (int x=0; x <j+i-taxa; ++x) cell_j=cell_j->next;
            }
            else {
                if (taxa==i+j+1) {
                    cell_j=present->distances;
                    continue;
                }
                else if (taxa==i+j+2);
                else cell_j=cell_j->next;
            }
            // add distance from new node to the taxa
            if ( taxa != i && taxa != i+j+1 ) {
                if (distant_taxa==0) distant_taxa = new_distance_node->distances;
                else {
                    distant_taxa->next = new distance_array;
                    distant_taxa = distant_taxa->next;
                }
                distant_taxa->next=0;
                distant_taxa->distance=(cell_i->distance + cell_j->distance - length)/2;
            }
            // remove cells for i and j
            if (taxa <i+j+1) {
                if (taxa < i) {
                    distance_array *prev;
                    if (cell_i == present->distances) prev = 0; // if we are in first column in row
                    else {
                        prev=present->distances;
                        while (prev->next != cell_i) prev=prev->next;
                    }
                    delete_distance_array_node(cell_i, prev, present);
                    cell_i=0;
                }
                if (taxa < j+i+1) {
                    distance_array *prev;
                    if (cell_j == present->distances) prev = 0;
                    else {
                        prev=present->distances;
                        while (prev->next != cell_j) prev=prev->next;
                    }
                    delete_distance_array_node(cell_j, prev, present);
                    cell_j=0;
                }
                present = present->next;
            }
        }
        node_and_distance_array* prev=distance_matrix;
        while (prev->next!=join_node1) prev=prev->next;
        delete_node_and_distance_array_node(join_node1, prev);
        prev=distance_matrix;
        while (prev->next!=join_node2) prev=prev->next;
        delete_node_and_distance_array_node(join_node2, prev);
        r = n_taxa_node_and_distance_array(distance_matrix);
    }
    root->left=distance_matrix->child;
    root->left->branchlength=0;
    root->left->parent=root;
    distance_matrix->child=0;
    root->right=distance_matrix->next->child;
    root->right->branchlength=distance_matrix->distances->distance;
    root->right->parent=root;
    distance_matrix->next->child=0;
    delete_node_and_distance_array(distance_matrix);
}

void njtree::delete_node_and_distance_array(node_and_distance_array* start_node) {
    if (start_node!=0) {
        if (start_node->next!=0) delete_node_and_distance_array(start_node->next);
        delete_distance_array(start_node->distances);
        start_node->distances=0;
        destroy_tree (start_node->child);
        start_node->child=0;
        delete start_node;
        if (start_node==distance_matrix) distance_matrix = 0;
    }
}

void njtree::delete_node_and_distance_array_node(node_and_distance_array* delete_node, node_and_distance_array* previous_node) {
    if (delete_node != 0) {
        delete_distance_array(delete_node->distances);
        if (previous_node == 0 && delete_node->next !=0) {
            destroy_tree (delete_node->child);
            delete_node->child=delete_node->next->child;
            delete_node->distances=delete_node->next->distances;
            delete_node->S=delete_node->next->S;
            delete_node_and_distance_array_node(delete_node->next,delete_node);
        }
        else {
            previous_node->next=delete_node->next;
            destroy_tree (delete_node->child);
            delete delete_node;
        }
    }
}
void njtree::delete_distance_array_node ( distance_array* delete_node, distance_array* previous_node, node_and_distance_array* father ) {
    if (delete_node != 0) {
        if (previous_node==0)  father->distances=delete_node->next;
        else previous_node->next = delete_node->next;
        delete delete_node;
    }
}

void njtree::delete_distance_array ( distance_array* start_node ) {
    if (start_node !=0) {
        if (start_node->next !=0) delete_distance_array( start_node->next );
        delete start_node->next;
    }
}

void njtree::read_distance_matrix( istream& infile, bool lables) {
    delete_node_and_distance_array( distance_matrix ); // clean up before new matrix
    distance_matrix = new node_and_distance_array; // prepare for new matrix
    distance_matrix->next=0;
    distance_matrix->child=0;
    distance_matrix->distances=0;
    distance_matrix->S=0;
    node_and_distance_array* present = distance_matrix;
    bool new_row = true;
    int n_taxa=0;
    string value;
    char temp;
    while (infile.good()) { // while we have input
        temp = infile.get(); // get input char by char
        /*if ((temp=='\n' || temp=='\r') && !new_row) { // if finding newline and not at start of row
            new_row=true;
        }*/
        if (temp == ' ' || temp=='\n' || temp=='\r' || temp=='\t' || infile.peek() == EOF) { // if finished reading value
            if (!value.empty()) { //if acctualy have value
                if (new_row) { // if at start of row
                    // if not in first row first collumn
                    if (present != distance_matrix || (present == distance_matrix && distance_matrix->child!=0 )) {
                        // create new row
                        present->next = new node_and_distance_array;
                        present->next->next=0;
                        present->next->S=0;
                        // go to new row
                        present=present->next;
                    }
                    present->child = new node; // add a node for the taxa of the row
                    present->child->left = 0;
                    present->child->right = 0;
                    present->child->parent = 0;
                    present->child->branchlength = 0;
                    if (lables) { // if lables are included
                        present->child->nodelabel = nodelabels.add_string(value);
                        present->distances = 0;
                    }
                    else { // if no lables first value is distance
                        stringstream converter;
                        converter << n_taxa;
                        present->child->nodelabel = nodelabels.add_string(converter.str()); // number for taxon name
                        present->distances = new distance_array;
                        present->distances->next = 0;
                        present->distances->distance = atof(value.c_str()); // add distance
                    }
                    ++n_taxa;
                    if (temp!='\n' && temp!='\r') new_row=false;
                }
                else if (present->distances == 0) { // if adding to first collumn in the row
                    present->distances = new distance_array;
                    present->distances->next = 0;
                    present->distances->distance = atof(value.c_str());
                }
                else add_to_distance_array (atof(value.c_str()), present->distances); // add distance
                value.clear();
            }
        }
        else {
            value += temp;
        }
        if (temp=='\n' || temp=='\r') new_row=true;
    }
    if (present->distances != 0) { // last taxa in the array should not have any distances since they are accounted for in above rows, so if last ataxa has distances add another taxa
        present->next = new node_and_distance_array;
        present->next->next=0;
        present->next->S=0;
        present=present->next;
        present->child = new node;
        present->child->left = 0;
        present->child->right = 0;
        present->child->parent = 0;
        present->child->branchlength = 0;
        stringstream converter;
        converter << n_taxa;
        present->child->nodelabel = nodelabels.add_string(converter.str());
        present->distances = 0;
    }
}

void njtree::print_node_and_distance_array ( node_and_distance_array* start_node, ostream& output ) {
    while ( start_node !=0 ) {
        if (start_node->child!=0) output << *start_node->child->nodelabel << ' ';
        if (start_node->distances!=0) print_distance_array( start_node->distances, output );
        output << endl;
        start_node = start_node->next;
    }
}

void njtree::print_distance_array ( distance_array* start_node, ostream& output ) {
    while ( start_node !=0 ) {
        output << start_node->distance << ' ';
        start_node = start_node->next;
    }
}

void njtree::add_to_distance_array ( float value, distance_array* start_node ) {
    if (start_node!=0) {
        while (start_node->next!=0) start_node=start_node->next;
        start_node->next = new distance_array;
        start_node->next->next=0;
        start_node->next->distance=value;
    }
}
/*
unsigned int njtree::n_nodes_in_array() {
    unsigned int n(0);
    node_and_distance_array* present = distance_matrix;
    while (present != 0) {
	present = present->next;
	n++;
    }
    return n;
}
*/
