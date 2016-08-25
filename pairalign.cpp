#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "seqpair.h"
#ifdef DATABASE
#include "sqlite3.h"
#endif //DATABASE
#include "align_group.h"
#ifdef PTHREAD
#include <pthread.h>
#endif /*PTHREAD*/
#include "seqdatabase.h"

using namespace std;

int n_threads=1;
#ifdef PTHREAD
pthread_mutex_t databasemutex;
#endif /*PTHREAD*/

//void pairalign ( istream* infile, const char output, const bool output_names, const bool aligned, const bool matrix );
void help();

// From alignment_groups

struct sequence_package {
    string accno1;
    string accno2;
    string sequence1;
    string sequence2;
    const string *table;
    seqdatabase* db;
    const float *cut_off;
    align_group *deviations;
    bool aligned;
    char output;
    bool output_names;
    bool matrix;
};

#ifdef PTHREAD
void *pthread_align_pair (void *threadarg);
#endif /* PTHREAD */

void align_pair ( sequence_package *two_sequences );
void cluster ( seqdatabase& db, const string table, const float cut_off, const int min_length, const bool only_lead, const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet );
void cluster_each_table ( const string& file, const char* databasetype, const string cut_off, const int min_length, const bool only_lead, const bool perform_clustering, const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet );

/********************************************************************
 **** MAIN FUNCTION, parse arguments and execute cluster function ****
********************************************************************/
int main (int argc, char *argv[]) {
    char output_mode = 'a';
    bool matrix = false;
    bool quiet(true);
    bool output_names = false;
    bool aligned = false;
    string file_name;
    ifstream infile;
    //istream* infile_stream = &std::cin;
    stringstream stdin_holder;
    // From alignment_groups
    time_t rawtime = time(0);
    struct tm * timeinfo = localtime( &rawtime );
    // Cut off value set to default
    string cut_off="all,0.999";
    int min_length = 100;
    bool only_lead = 0;
    bool perform_clustering = true;
    char db_type[7] = "fasta"; // alternative fasta
    //std::cout << "The program was called with the following command:" << endl;
    //for (int i=0; i<argc; ++i) std::cout << argv[i] << ' ';
    //std::cout << endl << endl;

   // Parse arguments
   for (int i = 1; i < argc; ++i) {
       if ( !strcmp(argv[i],"-a") || !strcmp(argv[i],"--alignments") ) output_mode = 'a';
       else if ( !strcmp(argv[i],"-d") || !strcmp(argv[i],"--distances") ) output_mode = 'd';
       else if ( !strcmp(argv[i],"-p") || !strcmp(argv[i],"--proportion_difference") ) output_mode = 'p';
       else if ( !strcmp(argv[i],"-s") || !strcmp(argv[i],"--similarity") ) output_mode = 's';
       else if ( !strcmp(argv[i],"-j") || !strcmp(argv[i],"--jc_distance") ) output_mode = 'j';
       else if ( !strcmp(argv[i],"-i") || !strcmp(argv[i],"--difference") ) output_mode = 'c';
       else if ( !strcmp(argv[i],"-m") || !strcmp(argv[i],"--matrix") ) {
	    matrix = true;
	    if (output_mode == 'a') output_mode = 'p';
	}
       else if ( !strcmp(argv[i],"-n") || !strcmp(argv[i],"--names") ) output_names = true;
       else if ( !strcmp(argv[i],"-A") || !strcmp(argv[i],"--aligned") ) aligned = true;
//       else if ( !strcmp(argv[i],"-q") || !strcmp(argv[i],"--quiet") ) quiet = true;
	else if (!strcmp(argv[i],"-v") || !strcmp(argv[i],"--verbose")) quiet = false;
       #ifdef PTHREAD
       else if ( !strcmp(argv[i],"-T") || !strcmp(argv[i],"--threads") ) {
           ++i;
           if (i < argc && argv[i][0] != '-') {
               n_threads = atoi(argv[i]);
               //std::cerr << "Number of threads set to: " << n_threads << endl;
           }
           else {
               std::cerr << "--threads or -T must be followed by a integer value, e.g. -T 4. Quiting quietly." << endl;
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
	else if ( !strcmp(argv[i],"-a") || !strcmp(argv[i],"--alignment_groups")) {
	    output_mode = 'A';
	}
	/// From alignment_groups
	else if ( !strcmp(argv[i],"-c") || !strcmp(argv[i],"--cut_off") ) {
	    ++i;
	    if (i < argc && argv[i][0] != '-') { 
		cut_off = argv[i];
		//std::cerr << "Cut offs for clustering set to: " << cut_off << endl;
	    }
	    else {
	       	std::cerr << "--cut_off or -c must be followed by a comma separated string. Quitting quietly." << endl;
		return 0; 
	    }
	}
	else if ( !strcmp(argv[i],"-m") || !strcmp(argv[i],"--min_length") ) {
	    ++i;
	    if (i < argc && argv[i][0] != '-') {
	     	min_length = atoi(argv[i]);
	      	//std::cerr << "Minimum sequence length to cluster set to: " << min_length << endl;
	    }
	    else {
	       	std::cerr << "--min_length or -m must be followed by a integer value, e.g. -m 100. Quitting quietly." << endl;
		return 0;
	    }
	}
	else if ( !strcmp(argv[i],"-P") || !strcmp(argv[i],"--previous_clusters") ) {
           only_lead = 1;
           //std::cerr << "Only lead sequences of previous clusters will be considered." << endl;
	}
	else if ( !strcmp(argv[i],"-O") || !strcmp(argv[i],"--no_cluster") ) {
	    perform_clustering = false;
	    //std::cerr << "No clustering will be performed. Only constructing alignmenmt groups." << endl;
	}
	#ifdef DATABASE
	else if ( !strcmp(argv[i],"-D") || !strcmp(argv[i],"--db_type") ) {
	    if (i+1 < argc && argv[i+1][0] != '-') {
		++i;
		if (!strcmp(argv[i], "sqlite")) strcpy(db_type, "sqlite");
		else if (!strcmp(argv[i], "fasta")) strcpy(db_type, "fasta");
		else {
		    cerr << argv[i] << " is not a valid option. -D/--db_type can only take the values sqlite or fasta." << endl;
		    return 0;
		}
	    }
	    else {
		cerr << "-D/--db_type require sqlite or fasta as next option depending on database type." << endl;
	    }
	}
	#endif //DATABASE
	////////////////
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
//    if (output_mode == 'A') { // alignmentgroups
	bool error_flag = 0; // flag to indicate if we should quit
	if ( min_length < 0 ) {
	    std::cerr << "Minimum sequence length to consider for clustering must be positive integer, e.g. 100." << endl;
	    error_flag = 1;
	}
	//#ifdef PTHREAD
	if (n_threads < 1) {
	    std::cerr << "Number of threads (--threads or -T) must be 1 or more (not " << n_threads << "), e.g. 4." << endl;
	    error_flag = 1;
	}
	//#endif /* PTHREAD */
	#ifdef DATABASE
	// If database can not be opened, give instructions then quit
	if (!quiet) std::cerr << "Using database: " << file_name << '.' << endl;
	sqlite3 *db;
	if ( sqlite3_open(file_name.c_str(), &db) != 0 ) {
	    std::cerr << "Could not open " << file_name << ". Last argument has to be a sqlite database." << endl;
	    error_flag = 1;
	}
	sqlite3_close(db);
	#endif //DATABASE
	// if any problems above, quit!
	if (error_flag) {
	    std::cerr << "Quitting quietly." << endl;
	    return 0;
	}

	if (!quiet) std::cerr << "Started at:" << endl << asctime( timeinfo );
	// Execute clustering, last argument is interpreted as the database file
	cluster_each_table ( file_name, db_type, cut_off, min_length, only_lead, perform_clustering, aligned, output_mode, output_names, matrix, quiet );
	// Calculate end time
	rawtime = time(0);
	timeinfo = localtime( &rawtime );
	if (!quiet) std::cerr << "Ended at:" << endl << asctime( timeinfo );
     	return 0; // return normally

/*    }
    else { // pairalign
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
	#endif // PTHREAD 
	// if any problems above, quit!
	if (error_flag) {
	    std::cout << "Quitting quietly." << endl;
	    return 0;
	}
	pairalign ( infile_stream, output_mode, output_names, aligned, matrix );
	return 0; // return normally
    }*/
}

/*** Function to print help ***/
void help() {
    std::cout << "This program will perform pairwise alignment of sequences given in fasta format through standard in." << endl; 
    std::cout << "alignmentgroups (c) Martin Ryberg" << endl << endl;
    std::cout << "Usage:" << endl << "pairalign [arguments] < inputfile.fasta" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--aligned / -A                  input file is aligned." << endl;
    std::cout << "--alignment_groups              this option will cluster sequences that are similar and find the most inclusive taxa in a" << endl; 
    std::cout << "                                    hierarchy that are alignable. It need the taxonomy given after a (the first) | in the" << endl;
    std::cout << "                                    sequence name. Each taxa in the hierarchy should be separated by a semicolon, with the" << endl;
    std::cout << "                                    highest rank first and then increasingly nested levels until the lowest known level for" << endl;
    std::cout << "                                    the sequence.";
    #ifdef DATABASE
    std::cout << " Alternatively it require a SQLite database file with one table named gb_data" << endl;
    std::cout << "                                    that have a column named accno that contain accession numbers for each sequence, and" << endl; 
    std::cout << "                                    one column named taxon_string that contain the taxonomic hierarchy for the sequence" << endl;
    std::cout << "                                    as given above.It should also have a table for each gene region to be tested. Each" << endl;
    std::cout << "                                    of these tables should have a column named accno with the accession for each sequence," << endl;
    std::cout << "                                    one column named sequence with the sequence, and one column named cluster (with empty" << endl;
    std::cout << "                                    as default value)" << endl;
    #else
    std::cout << endl;
    #endif //DATABASE
    std::cout << "--alignments / -a               output aligned sequences pairwise." << endl;
////////////////////////////
    std::cout << "--cut_off / -c [0-1]            sets the cut off value in pairwise similarity for clustering of each gene. Should be" << endl;
    std::cout << "                                    comma separated string with gene name first and value second. If all genes should" << endl;
    std::cout << "                                    have same cut off 'all' could be given as gene. E.g. -c all,0.99 or ITS,0.97,LSU,0.99." << endl;
    #ifdef DATABASE
    std::cout << "--db_type / -D [sqlite/fasta]   sets if the database is in sqlite or fasta format." << endl;
    #endif // DATABASE
/////////////////////////////////
    std::cout << "--difference / -i               output difference between the Jukes-Cantor (JC) distance and proportion different sites." << endl;
    std::cout << "--distances / -d                output proportion different sites, JC distance, and diference between the two." << endl;
    std::cout << "--help / -h                     print this help." << endl;
    std::cout << "--jc_distance / -j              output Jukes-Cantor (JC) distance." << endl;
    std::cout << "--matrix / -m                   output in the form of a space separated left-upper triangular matrix." << endl;
    std::cout << "--min_length / -m [integer]     sets the minimum length of sequences to consider for clustering, e.g. -m 100." << endl;
    std::cout << "--names / -n                    output sequence names (if outputing alignments then in fasta format)." << endl;
    std::cout << "--no_cluster / -O               turn clustering off when running alignment_groups. Only calculating alignment groups." << endl;
    #ifdef DATABASE
    std::cout << "--previous_clusters / -P        only sequences with cluster marked as 'lead' will be considered for further clustering." << endl;
    #endif // DATABASE
    std::cout << "--proportion_difference / -p    output proportion sites that are different." << endl;
//    std::cout << "--quiet / -q                    suppress additional output regarding the run of the program." << endl;
    std::cout << "--similarity / -s               output similarity between sequences (1-proportion different)." << endl;
    #ifdef PTHREAD
    std::cout << "--threads / -T [integer]        set the number of threads additional to the controling thread, e.g. -T 4." << endl;
    #endif /* PTHREAD */
    std::cout << "--verbose / -v                  get additional output." << endl;
}
/*
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
}*/
/// alignmentgroups
/*** Function to cluster each table ***/
void cluster_each_table ( const string& file, const char* databasetype, const string cut_off, const int min_length, const bool only_lead, const bool perform_clustering, const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet ) {
    seqdatabase database(file,databasetype);
    if (!database.alignment_groups_present())
	if (!database.create_alignment_groups()) return;
    vector<string> tables = database.tables_in_database();
    for (vector<string>::const_iterator table = tables.begin(); table != tables.end(); ++table) {
	if (!table->compare("gb_data") || !table->compare("alignments") || !table->compare("alignment_groups")) continue;
	float present_cut_off=0;
	if (perform_clustering) {
	    int length = cut_off.length();
	    string gene;
	    for (int i=0; i < length; ++i) {
		if (cut_off[i]==',') {
		    if(!gene.compare(*table) || !gene.compare("all")) {
			string number;
			++i;
			while (cut_off[i]!=',' && i<length) {
			    number += cut_off[i];
			    ++i;
			}
			present_cut_off = atof(number.c_str());
			break;
		    }
		}
		else gene+=cut_off[i];
	    }
	    if (present_cut_off < 0.000000001) {
		std::cerr << "Could not find appropriate cut off (" << present_cut_off<< ") for " << *table << ". Will only define alignment groups and not cluster." << endl;
		continue;
	    }
	}
	if (!quiet) std::cerr << "Checking " << *table << endl;
	cluster(database, *table, present_cut_off, min_length, only_lead, aligned, output, output_names, matrix, quiet);
	if (!strcmp(databasetype,"fasta")) {
	    database.print_clusters();
	}
    }
}

void cluster( seqdatabase& db, const string table, const float cut_off, const int min_length, const bool only_lead, const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet ) {
    #ifdef PTHREAD
    pthread_t thread[n_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    void *status;
    //string sequence1;
    //string accno1;
    //string cluster1;
    int next_thread=0;
    bool activated[n_threads];
    for (int i=0; i < n_threads; ++i) activated[i] = false;
    //#else /* PTHREAD */
    //sequence_package two_sequences;
    #endif /* PTHREAD */
    sequence_package two_sequences[n_threads];
    align_group deviations;
    //char mode='1';
    stringstream converter;
    converter << min_length;
    string char_min_length(converter.str());

    if (!quiet) std::cerr << "Starting pairwise alignment." << endl;
    string previous = "empty";
    //Need a fasta parser here as alternative, we should prepars
    //Maybe make an object that does the different steps and keep track
    unsigned int n_seq(0);
    if (db.initiate_sequence_retrieval(table, char_min_length)) {
	if (matrix && output == 'd') std::cout << "Proportion different/Similarity/Jukes-Cantor distance/Difference between JC and similarity" << endl; // from pairalign
    	while(!db.all_pairs()) {
	    bool new_row(false);
	    if (db.at_new_first())
		new_row = true;
    	    db.move_to_next_pair(only_lead);
	    if (matrix && new_row) {
		++n_seq;
	   	if (n_seq > 1) std::cout << endl;
       		if (output_names) std::cout << db.get_accno1() << ' ';
   		for (unsigned int i=1; i<n_seq; ++i) std::cout << ' ';
	    }
    	    #ifdef PTHREAD
	    if (n_threads > 1) {
		if (next_thread >= n_threads) next_thread = 0;
		int thread_code=0;
		if (activated[next_thread]) thread_code = pthread_join(thread[next_thread], &status);
		if (thread_code) {
		    std::cerr << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
		    exit (-1);
		}
		else {
		    two_sequences[next_thread].accno1 = db.get_accno1();
		    two_sequences[next_thread].sequence1 = db.get_sequence1();
		    two_sequences[next_thread].accno2 = db.get_accno2();
		    two_sequences[next_thread].sequence2 = db.get_sequence2();
		    two_sequences[next_thread].table = &table;
		    two_sequences[next_thread].db = &db;
		    two_sequences[next_thread].cut_off = &cut_off;
		    two_sequences[next_thread].deviations = &deviations;
		    two_sequences[next_thread].aligned = aligned;
		    two_sequences[next_thread].output = output;
		    two_sequences[next_thread].output_names = output_names;
		    two_sequences[next_thread].matrix = matrix;
		    thread_code = pthread_create(&thread[next_thread], &attr, pthread_align_pair, (void *) &two_sequences[next_thread]);
		    activated[next_thread] = true;
		    ++next_thread;
		    if (thread_code) {
			std::cerr << "ERROR!!! Return code from pthread_create() is: " << thread_code << endl;
			exit(-1);
		    }
		}
	    }
	    else {
    	    #endif /* PTHREAD */
    	    //#else /* PTHREAD */
    	    two_sequences[0].accno1 = db.get_accno1();
    	    two_sequences[0].sequence1 = db.get_sequence1();
    	    two_sequences[0].accno2 = db.get_accno2();
    	    two_sequences[0].sequence2 = db.get_sequence2();
    	    two_sequences[0].table = &table;
    	    two_sequences[0].db = &db;
    	    two_sequences[0].cut_off = &cut_off;
    	    two_sequences[0].deviations = &deviations;
	    two_sequences[0].aligned = aligned;
	    two_sequences[0].output = output;
	    two_sequences[0].output_names = output_names;
	    two_sequences[0].matrix = matrix;
    	    align_pair ( &two_sequences[0] ); // this is where the action is?
    	    //#endif /* PTHREAD */
    	    #ifdef PTHREAD
	    }
    	    #endif /* PTHREAD */
	    if (!quiet) cerr << '.';
	}
	if ( matrix && output_names ) std::cout << endl << db.get_accno2() << endl; // from pairalign
    }
    else cerr << "Could not initiate sequence retrieval. No aligning done for " << table << "." << endl;
    #ifdef PTHREAD
    pthread_attr_destroy(&attr);
    if (n_threads > 1) {
	for (int i=0; i < n_threads; ++i) {
	    int thread_code = pthread_join(thread[i], &status);
	    if (thread_code) {
		std::cerr << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
		exit (-1);
	    }
	}
    }
    #endif /* PTHREAD */
    if (!quiet) cerr << endl;
    if (output == 'A') {
	if (!quiet) std::cerr << "Finished aligning. Calculating mad to determine taxonomic level suitable for alignment." << endl;
	cout << "Alignment groups for " << table << endl;
	string alignment_groups = deviations.get_levels( table );
	std::cout << "    Alignment groups: " << alignment_groups << endl;
	//#if DATABASE
	db.insert_alignment_group(table, alignment_groups); 
	//#endif //DATABASE
	#ifdef DEBUG
	if (!quiet) cerr << "    Aprox. mad. entire group = " << deviations.aprox_mad() << endl;
	deviations.print_hierarchy();
	#endif /* DEBUG */
    }
}

#ifdef PTHREAD
void *pthread_align_pair (void *threadarg) {
    sequence_package *two_sequences = (sequence_package *) threadarg;
    align_pair (two_sequences);
    pthread_exit((void*)0);
}
#endif /* PROCESSTIME */

void align_pair ( sequence_package *two_sequences ) {
    seqpair sequences(two_sequences->sequence1,two_sequences->sequence2);
    #ifdef DEBUG
    cerr << two_sequences->accno1 << " " << two_sequences->accno2 << endl;
    #endif //DEBUG
    if (!two_sequences->aligned) {
	sequences.set_cost_matrix( 7, -5 );
	sequences.align();
    }
    #ifdef PTHREAD
    if (n_threads >1)
	pthread_mutex_lock (&databasemutex);
    #endif /* PTHREAD */
    if (two_sequences->output == 'A') {
	string taxon_string1 = two_sequences->db->get_taxon_string(two_sequences->accno1);
	string taxon_string2 = two_sequences->db->get_taxon_string(two_sequences->accno2);
	#ifdef DEBUG
	cerr << "Taxon string1: " << taxon_string1 << endl << "Taxon string2: " << taxon_string2 << endl;
	#endif //DEBUG
	if (*two_sequences->cut_off > 0.000000001) {
	    #ifdef DEBUG
	    cerr << "Staring to compare seq for clustering." << endl;
	    #endif //DEBUG
	    string cluster1 = two_sequences->db->get_cluster( two_sequences->accno1, *two_sequences->table );
	    string cluster2 = two_sequences->db->get_cluster( two_sequences->accno2, *two_sequences->table );
	    #ifdef DEBUG
	    cerr << "Seq1 belong to cluster: " << cluster1 << ". Seq2 belong to cluster: " << cluster2 << "." << endl;
	    #endif //DEBUG
	    if (sequences.similarity(true) > *two_sequences->cut_off) {
		#ifdef DEBUG
		cerr << "Clustering seq togather." << endl;
		#endif //DEBUG
		float comp1;
		float comp2;
		if (!cluster1.compare( "empty" ) || !cluster1.compare( "lead" )) comp1 = two_sequences->db->get_comp_value( two_sequences->accno1, *two_sequences->table );
		else comp1 = two_sequences->db->get_comp_value( cluster1, *two_sequences->table );
		if (!cluster2.compare( "empty" ) || !cluster2.compare( "lead" )) comp2 = two_sequences->db->get_comp_value( two_sequences->accno2, *two_sequences->table );
		else comp2 = two_sequences->db->get_comp_value( cluster2, *two_sequences->table );
		// If sequence 1 better
		if (comp1 >= comp2) {
		    #ifdef DEBUG
		    cerr << "Seq1 longer than seq2." << endl;
		    #endif //DEBUG
		    if ( !cluster1.compare( "empty" ) ) {
			two_sequences->db->clust_update( two_sequences->accno1, "lead", *two_sequences->table, 1 );
			two_sequences->db->clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, 1 );
			if (!cluster2.compare( "lead" ) || !cluster2.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, 0 );
			else {
			    two_sequences->db->clust_update( cluster2, two_sequences->accno1, *two_sequences->table, 1 );
			    two_sequences->db->clust_update( cluster2, two_sequences->accno1, *two_sequences->table, 0 );
			}
		    }
		    else if ( !cluster1.compare( "lead" ) ) {
			two_sequences->db->clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, 1 );
			if (!cluster2.compare( "lead" ) || !cluster2.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, 0 );
			else {
			    two_sequences->db->clust_update( cluster2, two_sequences->accno1, *two_sequences->table, 1 );
			    two_sequences->db->clust_update( cluster2, two_sequences->accno1, *two_sequences->table, 0 );
			}
		    }
		    else {
			two_sequences->db->clust_update( two_sequences->accno2, cluster1, *two_sequences->table, 1 );
			if (!cluster2.compare( "lead" ) || !cluster2.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno2, cluster1, *two_sequences->table, 0 );
			else {
			    two_sequences->db->clust_update( cluster2, cluster1, *two_sequences->table, 1 );
			    two_sequences->db->clust_update( cluster2, cluster1, *two_sequences->table, 0 );
			}
		    }
		}
		// If sequence 2 better
		else {
		    #ifdef DEBUG
		    cerr << "Seq2 longer than seq1." << endl;
		    #endif //DEBUG
		    if ( !cluster2.compare( "empty" ) ) { // if no previous annotation
			two_sequences->db->clust_update( two_sequences->accno2, "lead", *two_sequences->table, 1 ); // set best sequence to lead
			two_sequences->db->clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, 1 ); // set worse sequence to the same cluster
			if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, 0 ); // if worse sequence lead or empty uppdate sequences in that cluster
			else { // if not lead
			    two_sequences->db->clust_update( cluster1, two_sequences->accno2, *two_sequences->table, 1 ); // change lead of cluster to point to new best sequence
			    two_sequences->db->clust_update( cluster1, two_sequences->accno2, *two_sequences->table, 0 ); // change all sequences previusly pointing to that sequence to point to new best sequence
			}
		    }
		    else if ( !cluster2.compare( "lead" ) ) {
			two_sequences->db->clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, 1 ); // set worse sequence to the same cluster
			if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, 0 ); // if worse sequence lead or empty uppdate sequences in that cluster
			else {
			    two_sequences->db->clust_update( cluster1, two_sequences->accno2, *two_sequences->table, 1 ); // change lead of cluster to point to new best sequence
			    two_sequences->db->clust_update( cluster1, two_sequences->accno2, *two_sequences->table, 0 ); // change all sequences previusly pointing to that sequence to point to new best sequence
			}
		    }
		    else {
			two_sequences->db->clust_update( two_sequences->accno1, cluster2, *two_sequences->table, 1 ); // set worse sequence to point to the same sequence as better sequence
			if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno1, cluster2, *two_sequences->table, 0 ); // set worse sequence to point to same sequence as better sequence
			else {
			    two_sequences->db->clust_update( cluster1, cluster2, *two_sequences->table, 1 ); // change lead of worse cluster to point to the same sequence as new best sequence
			    two_sequences->db->clust_update( cluster1, cluster2, *two_sequences->table, 0 ); // change rest of cluster to point to the same sequence as new best sequence
			}
		    }
		}
		#ifdef DEBUG
		cerr << "Clustered seq togather." << endl;
		#endif //DEBUG
	    }
	    else {
		#ifdef DEBUG
		cerr << "Sequences are not clustered togather." << endl;
		#endif //DEBUG
		if ( !cluster1.compare( "empty" ) ) two_sequences->db->clust_update( two_sequences->accno1, "lead", *two_sequences->table, 1 );
		if ( !cluster2.compare( "empty" ) ) two_sequences->db->clust_update( two_sequences->accno2, "lead", *two_sequences->table, 1 );
		if ((!cluster1.compare( "empty" ) || !cluster1.compare( "lead" )) && (!cluster2.compare( "empty" ) || !cluster2.compare( "lead" ))) {
		    #ifdef DEBUG
		    cerr << "Since sequences may get included in alignment they are counted towards the MAD." << endl;
		    #endif
		    two_sequences->deviations->insert_value( taxon_string1, taxon_string2, sequences.jc_distance()-(1-sequences.similarity()) );
		}
	    }
	}
	else two_sequences->deviations->insert_value( taxon_string1, taxon_string2, sequences.jc_distance()-(1-sequences.similarity()) );
    }
    else {
	if ( two_sequences->output == 'a' ) {
	    if (two_sequences->output_names)
		std::cout << '>' << two_sequences->accno1 << endl;
	    std::cout << sequences.get_x() << endl;
	    if (two_sequences->output_names) 
		std::cout << '>' << two_sequences->accno2 << endl;
	    std::cout << sequences.get_y() << endl;
	}
	else if ( two_sequences->output == 'd' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) std::cout << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		std::cout << "Proportion sites that are different (and similarity): " << 1-sequences.similarity() << " (" << sequences.similarity() << "), Jukes-Cantor distance: " << sequences.jc_distance() << ", difference between the two: " << sequences.jc_distance()-(1.0-sequences.similarity()) << '.' << endl;
	    }
	    else {
		std::cout << 1-sequences.similarity() << "/" << sequences.similarity() << "/" << sequences.jc_distance() << "/" << sequences.jc_distance()-(1.0-sequences.similarity()) << '.';
	    }
	}
	else if ( two_sequences->output == 'c' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) std::cout << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		    std::cout << sequences.jc_distance()-(1.0-sequences.similarity()) << endl;
	    }
	    else std::cout << sequences.jc_distance()-(1.0-sequences.similarity()) << ' ';
	}
	else if ( two_sequences->output == 'j' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) std::cout << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		std::cout << sequences.jc_distance() << endl;
	    }
	    else std::cout << sequences.jc_distance() << ' ';
	}
	else if ( two_sequences->output == 'p' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) std::cout << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		std::cout << 1-sequences.similarity() << endl;
	    }
	    else std::cout << 1-sequences.similarity() << ' ';
	}
	else if ( two_sequences->output == 's' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) std::cout << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		std::cout << sequences.similarity() << endl;
	    }
	    else std::cout << sequences.similarity() << ' ';
	}
    }
    #ifdef PTHREAD
    if (n_threads > 1)
	pthread_mutex_unlock(&databasemutex);
    #endif /* PTHREAD */
}


