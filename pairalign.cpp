#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include "seqpair.h"
#ifdef PTHREAD
#include <pthread.h>
#endif /*PTHREAD*/

using namespace std;

#ifdef PTHREAD
int n_threads=4;
pthread_mutex_t databasemutex;
#endif /*PTHREAD*/

void pairalign ( istream* infile, const char output, const bool output_names, const bool aligned, const bool matrix );
void help();

/********************************************************************
 **** MAIN FUNCTION, parse arguments and execute cluster function ****
********************************************************************/
int main (int argc, char *argv[]) {
   char output_mode = 'a';
   bool matrix = false;
   bool quiet = false;
   bool output_names = false;
   bool aligned = false;
    string file_name;
    ifstream infile;
    istream* infile_stream = &std::cin;
    stringstream stdin_holder;
   // Parse arguments
   for (int i = 1; i < argc; ++i) {
       if ( !strcmp(argv[i],"-a") || !strcmp(argv[i],"--alignments") ) output_mode = 'a';
       else if ( !strcmp(argv[i],"-d") || !strcmp(argv[i],"--distances") ) output_mode = 'd';
       else if ( !strcmp(argv[i],"-p") || !strcmp(argv[i],"--proportion_difference") ) output_mode = 'p';
       else if ( !strcmp(argv[i],"-s") || !strcmp(argv[i],"--similarity") ) output_mode = 's';
       else if ( !strcmp(argv[i],"-j") || !strcmp(argv[i],"--jc_distance") ) output_mode = 'j';
       else if ( !strcmp(argv[i],"-c") || !strcmp(argv[i],"--difference") ) output_mode = 'c';
       else if ( !strcmp(argv[i],"-m") || !strcmp(argv[i],"--matrix") ) matrix = true;
       else if ( !strcmp(argv[i],"-n") || !strcmp(argv[i],"--names") ) output_names = true;
       else if ( !strcmp(argv[i],"-A") || !strcmp(argv[i],"--aligned") ) aligned = true;
       else if ( !strcmp(argv[i],"-q") || !strcmp(argv[i],"--quiet") ) quiet = true;
       #ifdef PTHREAD
       else if ( !strcmp(argv[i],"-T") || !strcmp(argv[i],"--threads") ) {
           ++i;
           if (i < argc && argv[i][0] != '-') {
               n_threads = atoi(argv[i]);
               std::cout << "Number of threads set to: " << n_threads << endl;
           }
           else {
               std::cout << "--threads or -T must be followed by a integer value, e.g. -T 4. Quiting quietly." << endl;
               return 0;
           }
       }
       #endif /* PTHREAD */
       else if ( !strcmp(argv[i],"-f") || !strcmp(argv[i],"--file") ) {
	    if (i+1 < argc && argv[i+1][0] != '-') {
		file_name = argv[++i];
	    }
	    else {
		cerr << "-f/--file needs to be followed by a file name." << endl;
		return 1;
	    }
	}
       else if ( !strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") ) {
           help();
           return 0;
       }
	else if ( i == argc-1 && argv[i][0] != '-' && file_name.empty() ) {
	    file_name = argv[i];
	}
       else if (i < argc) {
           std::cerr << "The program was called with the following command:" << endl;
           for (int j=0; j<argc; ++j) std::cerr << argv[j] << ' ';
           std::cerr << endl;
           std::cerr << "Argument " << argv[i] << " not recognized. For available arguments give -h or --help." << endl;
           return 0;
       }
   }
   if (!quiet) {
       std::cerr << "The program was called with the following command:" << endl;
       for (int i=0; i<argc; ++i) std::cerr << argv[i] << ' ';
       std::cerr << endl << endl;
   }
    if (!file_name.empty()) {
	infile.open(file_name.c_str(),std::ifstream::in);
	if (infile.good())
	    infile_stream = &infile;
	else {
	    cerr << "Could not open: " << file_name << endl;
	    return 1;
	}
    }
    else {
	while (*infile_stream) {
	    char temp;
	    *infile_stream >> noskipws >> temp;
	    stdin_holder<< temp;
	}
	infile_stream = &stdin_holder;
    }
   // Check if variables have reasonable values or else quit
   bool error_flag = 0; // flag to indicate if we should quit
   #ifdef PTHREAD
   if (n_threads < 1) {
       std::cout << "Number of threads (--threads or -T) must be more than 1, e.g. 4." << endl;
       error_flag = 1;
   }
   #endif /* PTHREAD */
   // if any problems above, quit!
   if (error_flag) {
       std::cout << "Quitting quietly." << endl;
       return 0;
   }
   pairalign ( infile_stream, output_mode, output_names, aligned, matrix );
   return 0; // return normally
}

/*** Function to print help ***/
void help() {
    std::cout << "This program will perform pairwise alignment of sequences given in fasta format through standard in." << endl; 
    std::cout << "alignmentgroups (c) Martin Ryberg" << endl << endl;
    std::cout << "Usage:" << endl << "pairalign [arguments] < inputfile.fasta" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--aligned / -A                  input file is aligned." << endl;
    std::cout << "--alignments / -a               output aligned sequences pairwise." << endl;
    std::cout << "--difference / -c               output difference between the Jukes-Cantor (JC) distance and proportion different sites." << endl;
    std::cout << "--distances / -d                output proportion different sites, JC distance, and diference between the two." << endl;
    std::cout << "--help / -h                     print this help." << endl;
    std::cout << "--jc_distance / -j              output Jukes-Cantor (JC) distance." << endl;
    std::cout << "--matrix / -m                   output in the form of a space separated left-upper triangular matrix." << endl;
    std::cout << "--names / -n                    output sequence names (if outputing alignments then in fasta format)." << endl;
    std::cout << "--proportion_difference / -p    output proportion sites that are different." << endl;
    std::cout << "--quiet / -q                    suppress additional output regarding the run of the program." << endl;
    std::cout << "--similarity / -s               output similarity between sequences (1-proportion different)." << endl;
    #ifdef PTHREAD
    std::cout << "--threads / -T [1+]             set the number of threads additional to the controling thread, e.g. -T 4." << endl;
    #endif /* PTHREAD */
}

void pairalign ( istream* infile, const char output, const bool output_names, const bool aligned, const bool matrix ) {
    string sequence1;
    string accno1;
    string sequence2;
    string accno2;
    char mode = 'n';
    int n_seq = 0;
    char inchar;
    int next_sequence;
    if (matrix && output == 'd') std::cout << "Proportion different/Similarity/Jukes-Cantor distance/Difference between JC and similarity" << endl;
    while (*infile) {
        *infile >> noskipws >> inchar;
        if (inchar == '>' || infile->peek()==EOF ) {
            if (inchar == '>' && (mode == 'a' || mode == 's')) {
                mode = 'b';
                ++n_seq;
                next_sequence = infile->tellg();
                if (matrix) {
                    if (n_seq > 1) std::cout << endl;
                    if (output_names) std::cout << accno1 << ' ';
                    for (int i=1; i<n_seq; ++i) std::cout << ' ';
                }
            }
            else if (mode == 'b' || mode == 't') {
		#ifdef DEBUG
		cerr << "Length seq1: " << sequence1.length() << "; seq2: " << sequence2.length() << endl;
		#endif //DEBUG
                seqpair sequences(sequence1,sequence2);
                if (!aligned) {
                    sequences.set_cost_matrix( 7, -5 );
                    sequences.align();
                }
                if ( output == 'a' ) {
                    if (output_names)
                        std::cout << '>' << accno1 << endl;
                    std::cout << sequences.get_x() << endl;
                    if (output_names) 
                        std::cout << '>' << accno2 << endl;
                    std::cout << sequences.get_y() << endl;
                }
                else if ( output == 'd' ) {
                    if (!matrix) {
                        if (output_names) std::cout << accno1 << " - " << accno2 << " | ";
                        std::cout << "Proportion sites that are different (and similarity): " << 1-sequences.similarity() << " (" << sequences.similarity() << "), Jukes-Cantor distance: " << sequences.jc_distance() << ", difference between the two: " << sequences.jc_distance()-(1.0-sequences.similarity()) << '.' << endl;
                    }
                    else {
                        std::cout << 1-sequences.similarity() << "/" << sequences.similarity() << "/" << sequences.jc_distance() << "/" << sequences.jc_distance()-(1.0-sequences.similarity()) << '.';
                    }
                }
                else if ( output == 'c' ) {
                    if (!matrix) {
                        if (output_names) std::cout << accno1 << " - " << accno2 << " | ";
                        std::cout << sequences.jc_distance()-(1.0-sequences.similarity()) << endl;
                    }
                    else std::cout << sequences.jc_distance()-(1.0-sequences.similarity()) << ' ';
                }
                else if ( output == 'j' ) {
                    if (!matrix) {
                        if (output_names) std::cout << accno1 << " - " << accno2 << " | ";
                        std::cout << sequences.jc_distance() << endl;
                    }
                    else std::cout << sequences.jc_distance() << ' ';
                }
                else if ( output == 'p' ) {
                    if (!matrix) {
                        if (output_names) std::cout << accno1 << " - " << accno2 << " | ";
                        std::cout << 1-sequences.similarity() << endl;
                    }
                    else std::cout << 1-sequences.similarity() << ' ';
                }
                else if ( output == 's' ) {
                    if (!matrix) {
                        if (output_names) std::cout << accno1 << " - " << accno2 << " | ";
                        std::cout << sequences.similarity() << endl;
                    }
                    else std::cout << sequences.similarity() << ' ';
                }
                sequence2.clear();
                accno2.clear();
                mode = 'b';
            }
            else if (mode == 'n') mode = 'a';
            if ( infile->peek()==EOF ) {
                if (matrix && output_names && (mode == 'a' || mode == 's')) std::cout << endl << accno1 << endl;
                else if (mode == 'b' || mode == 't') {
                    sequence1.clear();
                    accno1.clear();
                    mode = 'n';
                    infile->clear();
                    infile->seekg( next_sequence-1, ios::beg );
                }
            }
        }
        else if ( inchar == '\n' && mode == 'a' ) mode = 's';
        else if ( inchar == '\n' && mode == 'b' ) mode = 't';
        else if ( mode == 'a' ) accno1 += inchar;
        else if ( mode == 'b' ) accno2 += inchar;
        else if ( mode == 's' && (inchar == 'a' || inchar == 'A' || inchar == 'g' || inchar == 'G' 
                               || inchar == 'c' || inchar == 'C' || inchar == 't' || inchar == 'T' 
                               || inchar == 'r' || inchar == 'R' || inchar == 'y' || inchar == 'Y' 
                               || inchar == 's' || inchar == 'S' || inchar == 'w' || inchar == 'W' 
                               || inchar == 'k' || inchar == 'K' || inchar == 'm' || inchar == 'M' 
                               || inchar == 'b' || inchar == 'B' || inchar == 'd' || inchar == 'D' 
                               || inchar == 'h' || inchar == 'H' || inchar == 'v' || inchar == 'V' 
                               || inchar == 'n' || inchar == 'N' || inchar == '?' || (aligned && inchar == '-')) )
            sequence1 += inchar;
        else if ( mode == 't' && (inchar == 'a' || inchar == 'A' || inchar == 'g' || inchar == 'G'
                               || inchar == 'c' || inchar == 'C' || inchar == 't' || inchar == 'T'
                               || inchar == 'r' || inchar == 'R' || inchar == 'y' || inchar == 'Y'
                               || inchar == 's' || inchar == 'S' || inchar == 'w' || inchar == 'W'
                               || inchar == 'k' || inchar == 'K' || inchar == 'm' || inchar == 'M'
                               || inchar == 'b' || inchar == 'B' || inchar == 'd' || inchar == 'D'
                               || inchar == 'h' || inchar == 'H' || inchar == 'v' || inchar == 'V'
                               || inchar == 'n' || inchar == 'N' || inchar == '?' || (aligned && inchar == '-')) )
            sequence2 += inchar;
    }
}
