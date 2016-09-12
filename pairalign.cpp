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
    #if DEBUG
    cerr << "Debugging version!!!" << endl;
    #endif //DEBUG
    char output_mode = 'a';
    bool matrix = false;
    bool quiet(true);
    bool output_names = false;
    bool aligned = false;
    string file_name;
    ifstream infile;
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
	    }
	    else {
	       	std::cerr << "--min_length or -m must be followed by a integer value, e.g. -m 100. Quitting quietly." << endl;
		return 0;
	    }
	}
	else if ( !strcmp(argv[i],"-P") || !strcmp(argv[i],"--previous_clusters") ) {
           only_lead = 1;
	}
	else if ( !strcmp(argv[i],"-O") || !strcmp(argv[i],"--no_cluster") ) {
	    perform_clustering = false;
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
}

/*** Function to print help ***/
void help() {
    std::cout << "Pairalign " << VERSION << " will perform pairwise alignment of sequences given in fasta" << endl; 
    std::cout << "format through standard in." << endl;
    std::cout << "(c) Martin Ryberg " << YEAR << "." << endl << endl;
    std::cout << "Usage:" << endl << "pairalign [arguments] < inputfile.fasta" << endl << "pairalign [arguments] inputfile.fasta" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--aligned / -A                  input file is already aligned." << endl;
    std::cout << "--alignment_groups              this option will cluster sequences that are" << endl;
    std::cout << "                                similar and find the most inclusive taxa in a" << endl; 
    std::cout << "                                hierarchy that are alignable according to MAD" << endl;
    std::cout << "                                (Smith et al. 2009, BMC evol. Biol. 9:37). It" << endl;
    std::cout << "                                need the taxonomy given after a (the first) | in" << endl;
    std::cout << "                                the sequence name. Each taxa in the hierarchy" << endl;
    std::cout << "                                should be separated by a semicolon, with the" << endl;
    std::cout << "                                highest rank first and then increasingly nested" << endl;
    std::cout << "                                levels until the lowest known level for the" << endl;
    std::cout << "                                sequence.The groups that can be aligned are put" << endl;
    std::cout << "                                in a file with the ending .alignment_groups and" << endl;
    std::cout << "                                clusters in a file with ending .clusters." << endl; 
    #ifdef DATABASE
    std::cout << "                                Alternatively it require a SQLite database file" << endl;
    std::cout << "                                with one table named gb_data that have a column" << endl;
    std::cout << "                                named accno that contain accession numbers for" << endl; 
    std::cout << "                                each sequence, and one column named taxon_string" << endl;
    std::cout << "                                that contain the taxonomic hierarchy for the" << endl;
    std::cout << "                                sequence as given above. It should also have a" << endl;
    std::cout << "                                table for each gene region to be tested. Each of" << endl;
    std::cout << "                                these tables should have a column named accno" << endl;
    std::cout << "                                with the accession for each sequence, one column" << endl;
    std::cout << "                                named sequence with the sequence, and one column" << endl;
    std::cout << "                                named cluster (with empty as default value). The" << endl;
    std::cout << "                                cluster output will be saved to the database." << endl;
    #else
    std::cout << endl;
    #endif //DATABASE
    std::cout << "--alignments / -a               output aligned sequences pairwise." << endl;
////////////////////////////
    std::cout << "--cut_off / -c [value/s]        sets the cut off value in pairwise similarity" << endl;
    std::cout << "                                for clustering of each gene. Should be comma" << endl;
    std::cout << "                                separated string with gene name first and value" << endl;
    std::cout << "                                second. If all genes should have same cut off" << endl;
    std::cout << "                                'all' could be given as gene. E.g. -c all,0.99" << endl;
    std::cout << "                                or ITS,0.97,LSU,0.99." << endl;
    #ifdef DATABASE
    std::cout << "--db_type / -D [sqlite/fasta]   sets if the database is in sqlite or fasta" << endl;
    std::cout << "                                format." << endl;
    #endif // DATABASE
/////////////////////////////////
    std::cout << "--difference / -i               output difference between the Jukes-Cantor (JC)" << endl;
    std::cout << "                                distance and proportion different sites." << endl;
    std::cout << "--distances / -d                output proportion different sites, JC distance," << endl;
    std::cout << "                                and diference between the two." << endl;
    std::cout << "--help / -h                     print this help." << endl;
    std::cout << "--jc_distance / -j              output Jukes-Cantor (JC) distance." << endl;
    std::cout << "--matrix / -m                   output in the form of a space separated" << endl;
    std::cout << "                                left-upper triangular matrix." << endl;
    std::cout << "--min_length / -m [value]       sets the minimum length of sequences to consider"<< endl;
    std::cout << "                                for clustering, e.g. -m 100." << endl;
    std::cout << "--names / -n                    output sequence names (if outputting alignments" << endl;
    std::cout << "                                then in fasta format)." << endl;
    std::cout << "--no_cluster / -O               turn clustering off when running" << endl;
    std::cout << "                                --alignment_groups. Only calculating alignment" << endl;
    std::cout << "                                groups." << endl;
    #ifdef DATABASE
    std::cout << "--previous_clusters / -P        only consider sequences with cluster marked as" << endl;
    std::cout << "                                'lead' for further clustering." << endl;
    #endif // DATABASE
    std::cout << "--proportion_difference / -p    output proportion sites that are different." << endl;
//    std::cout << "--quiet / -q                    suppress additional output regarding the run of the program." << endl;
    std::cout << "--similarity / -s               output similarity between sequences (1-proportion" << endl;
    std::cout << "                                different)." << endl;
    #ifdef PTHREAD
    std::cout << "--threads / -T [number]         set the number of threads additional to the" << endl;
    std::cout << "                                controlling thread, e.g. -T 4. The order of" << endl;
    std::cout << "                                output between pairs are not guaranteed for more" << endl;
    std::cout << "                                than one thread (default 1)." << endl;
    #endif /* PTHREAD */
    std::cout << "--verbose / -v                  get additional output." << endl;
}
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
    #ifdef DEBUG
    cerr << "Starting my buisness" << endl;
    #endif //DEBUG

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
    vector<sequence_package> two_sequences;
    for (unsigned int i=0;i < n_threads; ++i) two_sequences.push_back(sequence_package());
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
    	    db.move_to_next_pair(only_lead);
	    if (db.at_new_first())
		new_row = true;
	    if (matrix && new_row) {
		++n_seq;
	   	if (n_seq > 1) std::cout << endl;
       		if (output_names) std::cout << db.get_accno1() << ' ';
   		for (unsigned int i=1; i<n_seq; ++i) std::cout << ' ';
	    }
    	    #ifdef PTHREAD
	    if (n_threads > 1) {
		#ifdef DEBUG
		cerr << "Prepairing pthread alignment." << endl;
		#endif //DEBUG
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
		#ifdef DEBUG
		cerr << "Prepairing alignment." << endl;
		#endif //DEBUG
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

