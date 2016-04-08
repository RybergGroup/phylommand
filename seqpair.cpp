#include "seqpair.h"

void seqpair::set_DNA_alphapet() {
    alphabet['A'].set(0);
    alphabet['G'].set(1);
    alphabet['C'].set(2);
    alphabet['T'].set(3);
    alphabet['-'].set(0); alphabet['-'].set(1); alphabet['-'].set(2); alphabet['-'].set(3);
    alphabet['R'].set(0); alphabet['R'].set(1);
    alphabet['Y'].set(2); alphabet['Y'].set(3);
    alphabet['S'].set(1); alphabet['S'].set(2);
    alphabet['W'].set(0); alphabet['W'].set(3);
    alphabet['K'].set(1); alphabet['K'].set(3);
    alphabet['M'].set(0); alphabet['M'].set(2);
    alphabet['B'].set(1); alphabet['B'].set(2); alphabet['B'].set(3);
    alphabet['D'].set(0); alphabet['D'].set(1); alphabet['D'].set(3);
    alphabet['H'].set(0); alphabet['H'].set(2); alphabet['H'].set(3);
    alphabet['V'].set(0); alphabet['V'].set(1); alphabet['V'].set(2);
    alphabet['N'].set(0); alphabet['N'].set(1); alphabet['N'].set(2); alphabet['N'].set(3);
    alphabet['.'].set(0); alphabet['.'].set(1); alphabet['.'].set(2); alphabet['.'].set(3);
}

string seqpair::translate_to_string ( const vector< bitset<SIZE> >& binary ) {
    string sequence;
    for (vector< bitset<SIZE> >::const_iterator i=binary.begin(); i != binary.end(); ++i) {
	for (map<char, bitset<SIZE> >::const_iterator j=alphabet.begin(); j!=alphabet.end(); ++j) {
	    if (*i == j->second) { sequence += j->first; break; }
	}
    }
    return sequence;
}

vector< bitset<SIZE> > seqpair::translate_to_binary ( const string& iupac ) {
    int n_bases = iupac.length();
    vector< bitset<SIZE> > sequence;
    for (int i=1; i<n_bases; ++i) {
	map<char, bitset<SIZE> >::const_iterator test = alphabet.find(iupac[i]);
	if (test != alphabet.end()) {
	    sequence.push_back(test->second);
	}
    }
    return sequence;
}


void seqpair::align( ) {
    const int n_bp_x=seq_x.size();
    const int n_bp_y=seq_y.size();
    int *aligned = new int [n_bp_x*n_bp_y];
    int *gap_y = new int [n_bp_x*n_bp_y];
    int *gap_x = new int [n_bp_x*n_bp_y];
    for (int i=0; i<n_bp_x; ++i) {
        for (int j=0; j<n_bp_y; ++j) {
            if (i==0 || j==0) {
                if (i==0 && j==0) {
                    aligned[i*n_bp_y+j] = cost(seq_x[i],seq_y[j]);
                    gap_y[i*n_bp_y+j] = GO;
                    gap_x[i*n_bp_y+j] = GO;
                }
                else if (i==0) {
                    aligned[i*n_bp_y+j] = gap_x[i*n_bp_y+j-1]+cost(seq_x[i],seq_y[j]);
                    gap_y[i*n_bp_y+j] = GO;
                    gap_x[i*n_bp_y+j] = gap_x[i*n_bp_y+j-1]+GE;
                }
                else if (j==0) {
                    aligned[i*n_bp_y+j] = gap_y[(i-1)*n_bp_y+j]+cost(seq_x[i],seq_y[j]);
                    gap_y[i*(n_bp_y)+j] = gap_y[(i-1)*n_bp_y+j]+GE;
                    gap_x[i*(n_bp_y)+j] = GO;

                }
            }
            else {
                if (aligned[(i-1)*n_bp_y+j-1] >= gap_y[(i-1)*n_bp_y+j] && aligned[(i-1)*n_bp_y+j-1] >= gap_x[i*n_bp_y+j-1] )  aligned[i*n_bp_y+j] = aligned[(i-1)*n_bp_y+j-1]+cost(seq_x[i],seq_y[j]);
                else if ( gap_y[(i-1)*n_bp_y+j] > aligned[(i-1)*n_bp_y+j-1] && gap_y[(i-1)*n_bp_y+j] > gap_x[i*n_bp_y+j-1] ) aligned[i*n_bp_y+j] = gap_y[(i-1)*n_bp_y+j]+cost(seq_x[i],seq_y[j]);
                else aligned[i*n_bp_y+j] = gap_x[i*n_bp_y+j-1]+cost(seq_x[i],seq_y[j]);
                if ( (aligned[(i-1)*n_bp_y+j-1]+GO) > (gap_y[(i-1)*n_bp_y+j]+GE) ) gap_y[i*n_bp_y+j] = aligned[(i-1)*n_bp_y+j-1]+GO;
                else gap_y[i*n_bp_y+j] = gap_y[(i-1)*n_bp_y+j]+GE;
                if ( (aligned[(i-1)*n_bp_y+j-1]+GO) > (gap_x[i*n_bp_y+j-1]+GE) ) gap_x[i*n_bp_y+j] = aligned[(i-1)*n_bp_y+j-1]+GO;
                else gap_x[i*n_bp_y+j] = gap_x[i*n_bp_y+j-1]+GE;
            }
        }
    }
    vector< bitset<SIZE> > seq_x_rev_aligned;
    vector< bitset<SIZE> > seq_y_rev_aligned;
    int i = n_bp_x-1;
    int j = n_bp_y-1;
    while (i>=0 || j>=0) {
        if ( i >=0 && j >=0 && aligned[i*n_bp_y+j] >= gap_y[i*n_bp_y+j] && aligned[i*n_bp_y+j] >= gap_x[i*n_bp_y+j] ) {
            seq_x_rev_aligned.push_back(seq_x[i]);
            seq_y_rev_aligned.push_back(seq_y[j]);
            --i;
            --j;
        }
        else if ( j < 0 || (i>=0 && gap_y[i*n_bp_y+j] >= aligned[i*n_bp_y+j] && gap_y[i*n_bp_y+j] >= gap_x[i*n_bp_y+j]) ) {
            bitset<SIZE> hyphen;
            seq_x_rev_aligned.push_back(seq_x[i]);
            seq_y_rev_aligned.push_back(hyphen);
            --i;
        }
        else if ( i < 0 || (j>=0 && gap_x[i*n_bp_y+j] >= aligned[i*n_bp_y+j] && gap_x[i*n_bp_y+j] >= gap_y[i*n_bp_y+j]) ) {
            bitset<SIZE> hyphen;
            seq_x_rev_aligned.push_back(hyphen);
            seq_y_rev_aligned.push_back(seq_y[j]);
            --j;
        }
    }

    delete [] aligned;
    delete [] gap_y;
    delete [] gap_x;
    i = seq_x_rev_aligned.size();
    seq_x.clear();
    while (i>0) seq_x.push_back(seq_x_rev_aligned[--i]);
    i = seq_y_rev_aligned.size(); 
    seq_y.clear();
    while (i>0) seq_y.push_back(seq_y_rev_aligned[--i]);
}

int seqpair::cost( const bitset<SIZE> bp_x, const bitset<SIZE> bp_y ) {
    int value = INT_MIN;
    for (unsigned int i=0; i < SIZE; ++i) {
        if (bp_x[i] || bp_y[i]) {
            for (int j=i; j < SIZE; ++j) {
                if ( (bp_x[i] && bp_y[j]) || (bp_y[i] && bp_x[j]) ) {
		    int coord = cost_matrix_coord(i,j);
                    if (cost_matrix[coord] > value) value = cost_matrix[coord];
                }
            }
        }
    }
    return value;
}

int seqpair::cost_matrix_coord(int x, int y) {
    if (y > x) {
	int t = y;
	y =x;
	x=t;
    }
    int c=0;
    for (int i=0; i<=x; ++i) c += i;
    c += y;
    return c;
}

void seqpair::set_cost_matrix ( int AA, int AG, int AC, int AT, int GG, int GC, int GT, int CC, int CT, int TT ) {
    cost_matrix [cost_matrix_coord(0,0)] = AA;
    cost_matrix [cost_matrix_coord(0,1)] = AG;
    cost_matrix [cost_matrix_coord(0,2)] = AC;
    cost_matrix [cost_matrix_coord(0,3)] = AT;
    cost_matrix [cost_matrix_coord(1,0)] = AG;
    cost_matrix [cost_matrix_coord(1,1)] = GG;
    cost_matrix [cost_matrix_coord(1,2)] = GC;
    cost_matrix [cost_matrix_coord(1,3)] = GT;
    cost_matrix [cost_matrix_coord(2,0)] = AC;
    cost_matrix [cost_matrix_coord(2,1)] = GC;
    cost_matrix [cost_matrix_coord(2,2)] = CC;
    cost_matrix [cost_matrix_coord(2,3)] = CT;
    cost_matrix [cost_matrix_coord(3,0)] = AT;
    cost_matrix [cost_matrix_coord(3,1)] = GT;
    cost_matrix [cost_matrix_coord(3,2)] = CT;
    cost_matrix [cost_matrix_coord(3,3)] = TT;
}

int seqpair::hamming_distance ( bool gap ) {
    int n = min (seq_x.size(),seq_y.size());
    int distance=0;
    for (int i=0; i < n; ++i) {
        if (!seq_x[i].any() && !seq_y[i].any()) continue;
        if (!gap && (!seq_x[i].any() || !seq_y[i].any())) continue;
	bitset<SIZE> test = seq_x[i] & seq_y[i];
        if (!test.any()) ++distance;
    }
    return distance;
}

double seqpair::similarity ( bool gap ) {
    int n_i = min (seq_x.size(),seq_y.size());
    int length = 0;
    for (int i=0; i<n_i; ++i) {
        if (seq_x[i]==0x00 && seq_y[i]==0x00) continue;
        else if (!gap && (seq_x[i]==0x00 || seq_y[i]==0x00)) continue;
        ++length;
    }
    if (length > 0) return 1.0-(hamming_distance(gap)/double(length));
    else return 1.0;
}

/*bool seqpair::same_iupac ( char bp_x, char bp_y ) {
    if (bp_x == 0x00 && bp_y == 0x00) return true;
    else if (bp_x & bp_y) return true;
}*/

