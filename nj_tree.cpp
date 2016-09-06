#include "nj_tree.h"

bool njtree::matrix_good() {
    unsigned int n=distance_matrix.size();
    for (vector<node_and_distance>::iterator present = distance_matrix.begin(); present != distance_matrix.end(); ++present) {
	if (present->distances.size() != n-1) return false;
	--n;
    }
    if (n != 0) return false;
    else return true;
}

void njtree::build_nj_tree (bool quiet) {
    int r = distance_matrix.size();//n_taxa_node_and_distance_array(distance_matrix);
    while (r>2) {
	#ifdef DEBUG
	print_node_and_distance_array(cerr);
        cerr << "Step 1: calculate sum of distances for each taxa" << endl;
	#endif //DEBUG
	for (vector<node_and_distance>::iterator present = distance_matrix.begin(); present != distance_matrix.end(); ++present)
	    present->S = 0;
	for (vector<node_and_distance>::iterator present = distance_matrix.begin(); present != distance_matrix.end(); ++present) {
	    vector<node_and_distance>::iterator pair = present;
	    for (vector<float>::iterator i=present->distances.begin(); i!=present->distances.end() && pair !=distance_matrix.end(); ++i) {
		present->S += *i;
		++pair;
		pair->S += *i;
	    }
	}
	#ifdef DEBUG
	for (vector<node_and_distance>::iterator present = distance_matrix.begin(); present != distance_matrix.end(); ++present)
	    cerr << "Sum: " << present->S << endl;
        cerr <<  "Step 2: calculate pair with smallest M" << endl;
	#endif //DEBUG
        float M=100000; // values in Q matrix
        int i=0;
        int j=0;
        int n=0;
	for (vector<node_and_distance>::iterator present = distance_matrix.begin(); present != distance_matrix.end(); ++present) {
            int n_taxa=0;
	    for (vector<float>::iterator I=present->distances.begin(); I!=present->distances.end(); ++I) {
                vector<node_and_distance>::iterator j_node = present+n_taxa+1;
                float value=(r-2) * *I - present->S - j_node->S;
		#ifdef DEBUG
		cerr << value << ' ';
		#endif //DEBUG
                if (value < M) { // No need to save all values just find lowest
                    M=value;
                    i=n;         // note the first
                    j=n_taxa;    // and second taxa in the pair
		    #ifdef DEBUG
		    cerr << "Join " << i << " and " << j << endl;
		    #endif //DEBUG
                }
                ++n_taxa;
            }
            ++n;
        }
	#ifdef DEBUG
	cerr << endl;
        cerr << "Step 3-4: Create a node joining selected taxa" << endl;
	#endif //DEBUG
        node* new_node = new node;
	#ifdef DEBUG
	cerr << "Created new node" << endl;
	#endif //DEBUG
        new_node->left = distance_matrix[i].child;
	#ifdef DEBUG
	cerr << "Linked taxa " << i << " (" << distance_matrix[i].child << ") to new node " << new_node;
       	if (new_node->left->nodelabel != 0) cerr << " (" << *new_node->left->nodelabel << ")";
       	cerr << endl; 
	#endif //DEBUG
        new_node->right = distance_matrix[i+j+1].child;
	#ifdef DEBUG
	cerr << "Linked taxa " << j+i+1 << " (" << distance_matrix[j+i+1].child << ") to new node " << new_node;
       	if (new_node->right->nodelabel != 0) cerr << " (" << *new_node->right->nodelabel << ")";
	cerr << endl; 
	#endif //DEBUG
        float length;
        length= distance_matrix[i].distances[j]; //cell->distance;
	#ifdef DEBUG
	cerr << "Distance between them: " << length << endl; 
	#endif //DEBUG
        new_node->left->branchlength = (length/2) + (distance_matrix[i].S - distance_matrix[j+i+1].S)/(2*(r-2));
	#ifdef DEBUG
	cerr << "Assigned branchlength to " << i << endl;
	#endif //DEBUG
        new_node->right->branchlength = length-new_node->left->branchlength;
	#ifdef DEBUG
	cerr << "Assigned branchlength to " << i+j+1 << endl;
	#endif //DEBUG
        new_node->left->parent = new_node;
        new_node->right->parent = new_node;
	#ifdef DEBUG
        cerr << "Step 5: calculate new distance matrix" << endl;
	#endif //DEBUG
        // create a new row
        node_and_distance new_distance_node;
        new_distance_node.child=new_node; // point it to new node
        new_distance_node.S=0;
        vector<float>::iterator cell_i;
        vector<float>::iterator cell_j;
        
        for (vector<node_and_distance>::iterator present = distance_matrix.begin();  present != distance_matrix.end(); ++present) { //int taxa=0; taxa<r && present != distance_matrix.end(); ++taxa) { // for each taxa
	    #ifdef DEBUG
	    cerr << "Working on taxon " << present - distance_matrix.begin() << endl;
	    #endif //DEBUG
            if (present - distance_matrix.begin() < i) { // if taxa before i the distance between them is in taxas row
                cell_i=present->distances.begin();
                for (int x=0; x <i-(present - distance_matrix.begin())-1 && cell_i != present->distances.end(); ++x)  ++cell_i;
            }
            else { // otherwise it is in is row
                if (present - distance_matrix.begin()==i) {
                    cell_i=present->distances.begin();
                    //++present;
                    continue;
                }
                else if (present - distance_matrix.begin()==i+1);
                else ++cell_i; // step one cell furter in the row
            }
	    #ifdef DEBUG
	    cerr << "Distance to first taxon: " << *cell_i << endl;
	    #endif //DEBUG
            if (present - distance_matrix.begin() < i+j+1) { // same for j
                cell_j=present->distances.begin();
                for (int x=0; x <j+i-(present - distance_matrix.begin()); ++x) ++cell_j;
            }
            else {
                if (present - distance_matrix.begin()==i+j+1) {
                    cell_j=present->distances.begin();
                    continue;
                }
                else if (present - distance_matrix.begin()==i+j+2);
                else ++cell_j;
            }
	    #ifdef DEBUG
	    cerr << "Distance to second taxon: " << *cell_j << endl;
	    #endif //DEBUG
            // add distance from new node to the taxa
            if ( present - distance_matrix.begin() != i && present - distance_matrix.begin() != i+j+1 ) {
                new_distance_node.distances.push_back((*cell_i + *cell_j - length)/2);
            }
            // remove cells for i and j
            if (present - distance_matrix.begin() < i+j+1) {
		bool i_removed(false);
                if (present - distance_matrix.begin() < i) {
		    present->distances.erase(cell_i);
                    cell_i=present->distances.end();
		    i_removed = true;
		    #ifdef DEBUG
		    cerr << "Deleted cell for first taxon, now " << present->distances.size() << " distances in row" << endl;
		    #endif //DEBUG
                }
                if (present - distance_matrix.begin() < j+i+1) {
		    if (i_removed) {
			cell_j=present->distances.begin();
			for (int x=0; x <j+i-(present - distance_matrix.begin())-1; ++x) ++cell_j;
		    }
		    present->distances.erase(cell_j);
                    cell_j=present->distances.end();
		    #ifdef DEBUG
		    cerr << "Deleted cell for second taxon, now " << present->distances.size() << " distances in row" << endl;
		    #endif //DEBUG
                }
            }
        }
	distance_matrix.erase(distance_matrix.begin()+i);
	distance_matrix.erase(distance_matrix.begin()+j+i);
	distance_matrix.insert(distance_matrix.begin(), new_distance_node);
        r = distance_matrix.size(); //n_taxa_node_and_distance_array(distance_matrix);
	if (!quiet) {
	    print_unfinished_tree(cerr);
	    cerr << endl;
	}
    }
    root->left=distance_matrix[0].child;
    root->left->branchlength=0;
    root->left->parent=root;
    distance_matrix[0].child=0;
    root->right=distance_matrix[1].child;
    root->right->branchlength=distance_matrix[0].distances[0];
    root->right->parent=root;
    distance_matrix[0].child=0;
    distance_matrix.clear();
    //delete_node_and_distance_array(distance_matrix);
}

/*void njtree::delete_node_and_distance_array(node_and_distance_array* start_node) {
    if (start_node!=0) {
        if (start_node->next!=0) delete_node_and_distance_array(start_node->next);
        //delete_distance_array(start_node->distances);
        //start_node->distances=0;
        destroy_tree (start_node->child);
        start_node->child=0;
        delete start_node;
        if (start_node==distance_matrix) distance_matrix = 0;
    }
}*/

/*void njtree::delete_node_and_distance_array_node(node_and_distance_array* delete_node, node_and_distance_array* previous_node) {
    if (delete_node != 0) {
        //delete_distance_array(delete_node->distances);
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
}*/
/*void njtree::delete_distance_array_node ( distance_array* delete_node, distance_array* previous_node, node_and_distance_array* father ) {
    if (delete_node != 0) {
        if (previous_node==0)  father->distances=delete_node->next;
        else previous_node->next = delete_node->next;
        delete delete_node;
    }
}*/

/*void njtree::delete_distance_array ( distance_array* start_node ) {
    if (start_node !=0) {
        if (start_node->next !=0) delete_distance_array( start_node->next );
        delete start_node->next;
    }
}*/


void njtree::read_distance_matrix( istream& infile, bool lables) {
    //delete_node_and_distance_array( distance_matrix ); // clean up before new matrix
    distance_matrix.clear();
    //distance_matrix = new node_and_distance_array; // prepare for new matrix
    //distance_matrix->next=0;
    //distance_matrix->child=0;
    //distance_matrix->distances=0;
    //distance_matrix->S=0;
    //node_and_distance_array* present = distance_matrix;
    bool new_row = true;
    int n_taxa=0;
    string value;
    char temp;
    while (infile.good()) { // while we have input
        temp = infile.get(); // get input char by char
        if (temp == ' ' || (temp=='\n' || temp=='\r') || temp=='\t' || infile.peek() == EOF) { // if finished reading value
            if (!value.empty()) { //if acctualy have value
                if (new_row) { // if at start of row
		    #ifdef DEBUG
		    cerr << "reading new row to the " << distance_matrix.size() << " that there already are" << endl;
		    #endif //DEBUG
                    // if not in first row first collumn
                    //if (present != distance_matrix || (present == distance_matrix && distance_matrix->child!=0 )) {
                    if (distance_matrix.empty() || distance_matrix.back().child!=0 ) {
			#ifdef DEBUG
			cerr << "Adding row to data matrix structure." << endl;
			#endif //DEBUG
                        // create new row
			distance_matrix.push_back(node_and_distance());
			//distance_matrix.back()->child = 0;
			distance_matrix.back().S = 0;
                        //present->next = new node_and_distance_array;
                        //present->next->next=0;
                        //present->next->S=0;
                        // go to new row
                        //present=present->next;
                    }
                    distance_matrix.back().child = new node; // add a node for the taxa of the row
                    distance_matrix.back().child->left = 0;
                    distance_matrix.back().child->right = 0;
                    distance_matrix.back().child->parent = 0;
                    distance_matrix.back().child->branchlength = 0;
                    if (lables) { // if lables are included
                        distance_matrix.back().child->nodelabel = nodelabels.add_string(value);
                        //present->distances = 0;
                    }
                    else { // if no lables first value is distance
                        stringstream converter;
                        converter << n_taxa;
                        distance_matrix.back().child->nodelabel = nodelabels.add_string(converter.str()); // number for taxon name
			distance_matrix.back().distances.push_back(atof(value.c_str()));
                        //present->distances = new distance_array;
                        //present->distances->next = 0;
                        //present->distances->distance = atof(value.c_str()); // add distance
                    }
                    ++n_taxa;
                    if (temp!='\n' && temp!='\r') new_row=false;
		    #ifdef DEBUGG
		    cerr << "Made new row" << endl;
		    #endif //DEBUGG
                }
		else { 
		    #ifdef DEBUG
		    cerr << "Adding distance " << atof(value.c_str()) << endl;
		    #endif //DEBUG
		    distance_matrix.back().distances.push_back(atof(value.c_str()));
	       	}
                /*else if (present->distances == 0) { // if adding to first collumn in the row
                    present->distances = new distance_array;
                    present->distances->next = 0;
                    present->distances->distance = atof(value.c_str());
                }*/
                //else add_to_distance_array (atof(value.c_str()), present->distances); // add distance
                value.clear();
            }
        }
        else {
            value += temp;
        }
        if (temp=='\n' || temp=='\r') new_row=true;
    }
    if (!distance_matrix.back().distances.empty()) { // != 0) { // all distances should already have been given for last taxa
	#ifdef DEBUG
	cerr << "Adding a last taxa" << endl;
	#endif //DEBUG
	distance_matrix.push_back(node_and_distance());
        //present->next = new node_and_distance_array;
        //present->next->next=0;
        distance_matrix.back().S=0;
        //present=present->next;
        distance_matrix.back().child = new node;
        distance_matrix.back().child->left = 0;
        distance_matrix.back().child->right = 0;
        distance_matrix.back().child->parent = 0;
        distance_matrix.back().child->branchlength = 0;
        stringstream converter;
        converter << n_taxa;
        distance_matrix.back().child->nodelabel = nodelabels.add_string(converter.str());
        //present->distances = 0;
    }
}

void njtree::print_node_and_distance_array ( vector<node_and_distance>::const_iterator start_node, ostream& output ) {
    //while ( start_node !=0 ) {
    for (;start_node != distance_matrix.end(); ++start_node) {
        if (start_node->child!=0 && start_node->child->nodelabel != 0) output << *start_node->child->nodelabel << ' ';
	else output << ' ';
        //if (start_node->distances!=0) print_distance_array( start_node->distances, output );
	for (vector<float>::const_iterator i = start_node->distances.begin(); i != start_node->distances.end(); ++i) output << *i << ' ';
        output << endl;
        //start_node = start_node->next;
    }
}

/*void njtree::print_distance_array ( distance_array* start_node, ostream& output ) {
    while ( start_node !=0 ) {
        output << start_node->distance << ' ';
        start_node = start_node->next;
    }
}*/

/*void njtree::add_to_distance_array ( float value, distance_array* start_node ) {
    if (start_node!=0) {
        while (start_node->next!=0) start_node=start_node->next;
        start_node->next = new distance_array;
        start_node->next->next=0;
        start_node->next->distance=value;
    }
}*/

/*unsigned int njtree::n_nodes_in_array() {
    unsigned int n(0);
    node_and_distance_array* present = distance_matrix;
    while (present != 0) {
	present = present->next;
	n++;
    }
    return n;
}*/

void njtree::print_unfinished_tree( ostream& output ) {
    output << '(';
    for (vector<node_and_distance>::const_iterator i=distance_matrix.begin(); i != distance_matrix.end(); ++i) {
	if (i != distance_matrix.begin()) output << ',';
	print_newick_subtree(output, i->child, -1, true, false );
    }
    cout << ");" << endl;
}
