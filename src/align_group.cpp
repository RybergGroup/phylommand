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

#include "align_group.h"
/*** This function find the most specific inclusive taxon for two sequences and call ****
**** find_node_insert_value to insert the value at that level in the hierarchy        ***/
void align_group::insert_value ( const string taxon_string1, const string taxon_string2, float value ) {
    #ifdef DEBUG
    cerr << "Begining to insert value" << endl << taxon_string1 << endl << taxon_string2 << endl;
    #endif //DEBUG
    string inclusive_taxon;
    int i=0;
    int j=0;
    int length1 = taxon_string1.length();
    int length2 = taxon_string2.length();
    #ifdef DEBUG
    cerr << "Length string1: " << length1 << " length string2: " << length2 << endl;
    #endif //DEBUG
    while (i < length1 && j < length2) {
	#ifdef DEBUG
	cerr << "Pos string1: " << i << " pos string2: " << j << endl;
	#endif //DEBUG
        string taxon1;
        string taxon2;
        while (taxon_string1[i] != ';' && i < length1) {
            taxon1 += taxon_string1[i++];
        }
        i += 2;
        while (taxon_string2[j] != ';' && j < length2) {
            taxon2 += taxon_string2[j++];
        }
        j += 2;
	#ifdef DEBUG
	cerr << "Taxon1: " << taxon1 << ", taxon2: " << taxon2<< endl;
	#endif //DEBUG
        if (!taxon1.compare(taxon2)) {
            inclusive_taxon += taxon1;
            inclusive_taxon += ';';
        }
        else { break; }
    }
    if ( inclusive_taxon.empty() ) inclusive_taxon = root->taxon;
    find_node_insert_value ( inclusive_taxon, value, root );
    #ifdef DEBUG
    cerr << "Finished insert value" << endl;
    #endif //DEBUG
}

/*** This node will find the right place in the hierarchy to insert a value, starting ****
**** from the root and given a taxon. If the taxon is not already in the hierarchy it ****
**** will be inserted.                                                               ***/
void align_group::find_node_insert_value ( string taxon, float value, node *leaf ) {
    #ifdef DEBUG
    cerr << "Finding and inserting value" << endl;
    #endif //DEBUG
    string heighest_taxon;
    string next_taxon;
    string rest;
    unsigned int i(0);
    while (taxon[i] != ';' && i < taxon.length()) heighest_taxon += taxon[i++];
    ++i;
    while (taxon[i] != ';' && i < taxon.length()) next_taxon += taxon[i++];
    while (i < taxon.length()) rest += taxon[i++];
    if (leaf->taxon.empty()) leaf->taxon.assign( heighest_taxon );
    if (!heighest_taxon.compare(leaf->taxon)) {
        if (next_taxon.empty()) {
            value *= precision;
            if (value > array_length-1) value=array_length-1;
            leaf->saturation_values[int(value)]+=1;
        }
        else {
            if ( leaf->nodes == 0 ) {
                leaf->nodes = new list<node>;
                add_node ( leaf->nodes, leaf, next_taxon );
                find_node_insert_value( next_taxon.append(rest), value, leaf->nodes->daughter );
            }
            else {
                list<node> *present = leaf->nodes;
                bool flag=0;
                while (1) {
                    if (!next_taxon.compare(present->daughter->taxon)) {
                        find_node_insert_value( next_taxon.append(rest), value, present->daughter );
                        flag = 1;
                    }
                    if (present->next == 0) break;
                    else present = present->next;
                }
                if (!flag) {
                    present->next = new list<node>;
                    add_node ( present->next, leaf, next_taxon );
                    find_node_insert_value( next_taxon.append(rest), value, present->next->daughter );
                }
            }
        }
	#ifdef DEBUG
	cerr << "Found and inserted value" << endl;
	#endif //DEBUG
    }
    else {
        std::cerr << "WARNING!!! Error in align_group::insert_value!!! Reached unexpected node!!! Value not inserted!!!" << endl;
    }
}

/*** This function will add a node/taxon to the hierarchical tree ***/
void align_group::add_node ( list<node> *position, node *parent, string taxon ) { 
    position->next = 0;
    position->daughter = new node; 
    position->daughter->nodes = 0;
    position->daughter->parent = parent;
    position->daughter->taxon = taxon;
    for (int i=0; i<array_length; ++i) position->daughter->saturation_values[i]=0;
}

/*** This function will go through the hierarchy from the root and until   ****
**** it hit a level that is alignable. It returns a semicolon separated  ****
**** string with the most specific taxa that are alignable based on MAD. ***/
string align_group::get_levels ( node *leaf, string gene ) { //sqlite3 *db, string gene ) {
    int values[array_length];
    for (int i=0; i<array_length; ++i) values[i]=0;
    add_values( values, leaf ); // add up the values of lower taxa
    #ifdef DEBUG
    std::cout << "Value array for " << leaf->taxon << ":" << endl;
    for (int i=0; i<array_length; ++i) std::cout << values[i] << ' ';
    std::cout << endl;
    #endif /* DEBUG */
    float mad = calc_aprox_mad( values );
    string taxon_string;
    if (mad < 0.01) {
        taxon_string=leaf->taxon;
        if (taxon_string.empty()) {
            return "empty";
        }
        return taxon_string.append("_A;");
    }
    // Otherwise go to lower taxa
    else if (leaf->nodes != 0) {
        list<node> *present=leaf->nodes;
        while (present!=0) { 
            // Add all the alignment groups from lower taxa
            taxon_string += get_levels(present->daughter, gene); //db, gene);
            present = present->next;
        }
        return taxon_string;
    }
    // If at tip return taxa that is not alignable mark it with T
    else {
        taxon_string=leaf->taxon;
        return taxon_string.append("_T;");
    }
}

/*** Returns the number of sequences for a specific taxon. Not used ****
**** any more.                                                      ***
int align_group::get_n_taxa( string taxon, sqlite3 *db, string table ) {
    sqlite3_stmt *statement;
    string query = "SELECT COUNT(accno) FROM ";
    query += table;
    query += " INNER JOIN ";
    query += gb_data;
    query += " ON ";
    query += table;
    query += ".accno = ";
    query += gb_data;
    query += ".accno WHERE ";
    query += gb_data;
    query += ".taxon_string='%";
    query += taxon;
    query += "%';"; // WARNING!!! This may include sequences from more inclusive taxa if they are based on the taxon name
    int n_taxa=0;
    if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
        if (sqlite3_step(statement) == SQLITE_ROW) {
            n_taxa=sqlite3_column_int(statement,0);
        }
    }
    sqlite3_finalize(statement);
    return n_taxa;
}
*/
/*** This function will calculate the mad value based on the value string ***/
float align_group::calc_aprox_mad( int values[] ) {
    int sum=0;
    int median;
    int deviation[array_length];
    int target;
    for (int i=0; i<array_length; ++i) deviation[i]=0;
    for (int i=0; i<array_length; ++i) sum += values[i];
    target = (sum/2)+sum%2;
    for (int i=0; i<array_length; ++i) {
        median = i;
        target -= values[i];
        if (target <= 0) break;
    }
    for (int i=0; i<array_length; ++i) deviation[abs(i-median)]+=values[i];
        target = (sum/2)+sum%2;
    for (int i=0; i<array_length; ++i) {
        median = i;
        target -= deviation[i];
        if (target <= 0) break;
    }
    return 1.4826*float(median)/precision;
}

/*** This function will sum the values for each entry ****
**** from taxa lower in the hierarchy.                 ***/
void align_group::add_values ( int values[], node *leaf ) {
    list<node> *present=leaf->nodes;
    while (present!=0) {
        add_values(values, present->daughter);
        present = present->next;
    }
    for (int i=0; i<array_length; ++i) values[i]+=leaf->saturation_values[i];
}

void align_group::print_hierarchy ( node *leaf ) {
    std::cout << leaf->taxon << endl;
    for (int i=0; i<array_length; ++i) std::cout << leaf->saturation_values[i] << ' ';
    std::cout << endl;
    list<node> *present=leaf->nodes;
    while (present!=0) {
        print_hierarchy (present->daughter);
        present = present->next;
    }
}
