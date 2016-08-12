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

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <bitset>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include <climits>

using namespace std;

template<int i> class size_triang_matrix {
    public:
    static const int result = i+size_triang_matrix<i-1>::result;
};
template<> class size_triang_matrix<1> {
    public:
    static const int result = 1;
};

/*** Class for phylogenetic trees ***/
class seqpair {
    public:
    seqpair () {
        GO = -10;
        GE = -1;
        set_cost_matrix();
	set_DNA_alphapet();

    };
    seqpair (string x, string y) {
	set_DNA_alphapet();
        translate_to_binary ( x, seq_x );
        translate_to_binary ( y, seq_y );
        GO = -10;
        GE = -1;
        set_cost_matrix();
    };
    void align();
    void add_bp_x ( const string& x ) { add_bp(seq_x, x); };
    void add_bp_y ( const string& y ) { add_bp(seq_y, y); };
    void add_chunk_x ( const string& x ) { add_bp(seq_x, x); };
    void add_chunk_y ( const string& y ) { add_bp(seq_y, y); };
    void set_GO ( int x ) { GO = x; };
    void set_GE ( int x ) { GE = x; };
    void set_cost_matrix ( int AA, int AG, int AC, int AT, int GG, int GC, int GT, int CC, int CT, int TT );
    void set_cost_matrix ( ) {
        set_cost_matrix ( 10, -1, -3, -4, 7, -5, -3, 9, 0, 8 );
    };
    void set_cost_matrix ( int match, int missmatch ) {
        set_cost_matrix ( match, missmatch, missmatch, missmatch, match, missmatch, missmatch, match, missmatch, match );
    };
    void set_gap_penalty ( int O, int E ) {
        GO = O;
        GE = E;
    };
    void set_DNA_alphapet();
    void set_alphapet( map<char, vector<int> >& translate );
    void set_alphapet( map<char, bitset<SIZE> >& translate );
    void print () { cout << translate_to_string (seq_x) << endl << translate_to_string (seq_y); };
    string get_x() { return translate_to_string ( seq_x ); };
    string get_y() { return translate_to_string ( seq_y ); };
    // hamming distance returns distance in number of bp
    int hamming_distance ( bool gap ); //if gapped is FALSE gaped sites are excluded from comparison
    int hamming_distance ( ) {
        return hamming_distance(0);
    };
    // similarity returns the proportion of identical sites
    double similarity ( bool gap ); //if gapped is FALSE gaped sites are excluded from comparison
    double similarity ( ) {
        return similarity(0);
    };
    // jc distance returns the jukes-cantor distance
    double jc_distance ( ) {
        double p=1-similarity(0);
        return log(1.0-(4.0/3.0)*p)*(-3.0/4.0);
    };

    private:
    vector< bitset<SIZE> > seq_x;
    vector< bitset<SIZE> > seq_y;
    string translate_to_string ( const vector< bitset<SIZE> >& binary );
    void translate_to_binary ( const string& iupac, vector< bitset<SIZE> >& sequence );
    map <char, bitset<SIZE> > alphabet;
    int GO;
    int GE;
    int cost_matrix_coord(int x, int y);
    int cost_matrix[size_triang_matrix<SIZE>::result]; // left lower triangle
    //map< pair<bitset<SIZE>, bitset<SIZE> >, int> cost_matrix;
    void add_bp ( vector< bitset<SIZE> >& seq, const string& x ) {
	vector< bitset<SIZE> > add;
	translate_to_binary ( x, add );
	for (vector< bitset<SIZE> >::iterator i=add.begin(); i != add.end(); ++i)
	    seq.push_back(*i);
    };
    int cost( const bitset<SIZE> bp_x, const bitset<SIZE> bp_y );
    //bool same_iupac ( char bp_x, char bp_y );
};
