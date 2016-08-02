#include <iostream>
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
#ifdef PTHREAD // if compiling multithreads
#include <pthread.h>
#endif /*PTHREAD*/
#include "seqdatabase.h"

using namespace std;

#ifdef PTHREAD
int n_threads=4;
pthread_mutex_t databasemutex;
#endif /*PTHREAD*/

struct sequence_package {
    string accno1;
    string accno2;
    string sequence1;
    string sequence2;
    const string *table;
    seqdatabase* db;
    const float *cut_off;
    align_group *deviations;
};

#ifdef PTHREAD
void *pthread_align_pair (void *threadarg);
#endif /* PTHREAD */

void align_pair ( sequence_package *two_sequences );
void cluster ( seqdatabase& db, const string table, const float cut_off, const int min_length, bool only_lead );
void cluster_each_table ( const char* file, const char* databasetype, const string cut_off, const int min_length, bool only_lead, bool perform_clustering );
void help();

/********************************************************************
**** MAIN FUNCTION, parse arguments and execute cluster function ****
********************************************************************/
int main (int argc, char *argv[]) {
    // Variables to calculate starting and end time
    time_t rawtime = time(0);
    struct tm * timeinfo = localtime( &rawtime );
    // Cut off value set to default
    string cut_off="all,0.999";
    int min_length = 100;
    bool only_lead = 0;
    bool perform_clustering = true;
    char db_type[7] = "fasta"; // alternative fasta
    std::cout << "The program was called with the following command:" << endl;
    for (int i=0; i<argc; ++i) std::cout << argv[i] << ' ';
    std::cout << endl << endl;
    // Parse arguments
    for (int i = 1; i < argc; ++i) {
	if ( !strcmp(argv[i],"-c") || !strcmp(argv[i],"--cut_off") ) {
	    ++i;
	    if (i < argc && argv[i][0] != '-') { 
		cut_off = argv[i];
		std::cout << "Cut offs for clustering set to: " << cut_off << endl;
	    }
	    else {
	       	std::cout << "--cut_off or -c must be followed by a comma separated string. Quitting quietly." << endl;
		return 0; 
	    }
	}
	else if ( !strcmp(argv[i],"-m") || !strcmp(argv[i],"--min_length") ) {
	    ++i;
	    if (i < argc && argv[i][0] != '-') {
	     	min_length = atoi(argv[i]);
	      	std::cout << "Minimum sequence length to cluster set to: " << min_length << endl;
	    }
	    else {
	       	std::cout << "--min_length or -m must be followed by a integer value, e.g. -m 100. Quitting quietly." << endl;
		return 0;
	    }
	}
	else if ( !strcmp(argv[i],"-p") || !strcmp(argv[i],"--previous_clusters") ) {
           only_lead = 1;
           std::cout << "Only lead sequences of previous clusters will be considered." << endl;
	}
	else if ( !strcmp(argv[i],"-n") || !strcmp(argv[i],"--no_cluster") ) {
	    perform_clustering = false;
	    std::cout << "No clustering will be performed. Only constructing alignmenmt groups." << endl;
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
	#endif // DATABASE
	#ifdef PTHREAD
       else if ( !strcmp(argv[i],"-T") || !strcmp(argv[i],"--threads") ) {
           if (i < argc && argv[i+1][0] != '-') {
               n_threads = atoi(argv[++i]);
               std::cout << "Number of threads set to: " << n_threads << endl;
           }
           else {
               std::cout << "--threads or -T must be followed by a integer value, e.g. -T 4. Quitting quietly." << endl;
               return 0;
           }
       }
       #endif /* PTHREAD */
       else if ( !strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") ) {
           help();
           return 0;
       }
       else if (i < argc-1) {
           std::cout << "Argument " << argv[i] << " not recognized. For available arguments give -h or --help." << endl;
           return 0;
       }
   }
   // Check if variables have reasonable values or else quit
   bool error_flag = 0; // flag to indicate if we should quit
   if ( min_length < 0 ) {
       std::cout << "Minimum sequence length to consider for clustering must be positive integer, e.g. 100." << endl;
       error_flag = 1;
   }
   #ifdef PTHREAD
   if (n_threads < 1) {
       std::cout << "Number of threads (--threads or -T) must be more than 1, e.g. 4." << endl;
       error_flag = 1;
   }
   #endif /* PTHREAD */
    #ifdef DATABASE
   // If database can not be opened, give instructions then quit
   std::cout << "Using database: " << argv[argc-1] << '.' << endl;
   sqlite3 *db;
   if ( sqlite3_open(argv[argc-1], &db) != 0 ) {
       std::cout << "Could not open " << argv[argc-1] << ". Last argument has to be a sqlite database." << endl;
       error_flag = 1;
   }
   sqlite3_close(db);
    #endif //DATABASE
   // if any problems above, quit!
   if (error_flag) {
       std::cout << "Quitting quietly." << endl;
       return 0;
   }

   std::cout << "Started at:" << endl << asctime( timeinfo );
   // Execute clustering, last argument is interpreted as the database file
   cluster_each_table ( argv[argc-1], db_type, cut_off, min_length, only_lead, perform_clustering );
   // Calculate end time
   rawtime = time(0);
   timeinfo = localtime( &rawtime );
   std::cout << "Ended at:" << endl << asctime( timeinfo );
   return 0; // return normally
}

/*** Function to print help ***/
void help() {
    std::cout << "This program will cluster sequences that are similar and and find the most inclusive taxa in a" << endl; 
    std::cout << "hierarchy that are alignable. It require a SQLite database file with one table named gb_data" << endl;
    std::cout << "that have a column named accno that contain accession numbers for each sequence, and one column" << endl; 
    std::cout << "named taxon_string that contain all taxa in the hierarchy that the sequence belong to separated" << endl;
    std::cout << "by semicolon. It should also have a table for each gene region to be tested. Each of these tables" << endl;
    std::cout << "should have a column named accno with the accession for each sequence, one column named sequence" << endl;
    std::cout << "with the sequence, and one column named cluster (with empty as default value)" << endl << endl;
    std::cout << "alignmentgroups version 0.3 (c) Martin Ryberg" << endl << endl;
    std::cout << "Usage:" << endl << "alignmentgroups [arguments] databasefile" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--cut_off / -c [0-1]          Sets the cut off value in pairwise similarity for clustering of each gene. Should be" << endl;
    std::cout << "                                  comma separated string with gene name first and value second. If all genes should" << endl;
    std::cout << "                                  have same cut off 'all' could be given as gene. E.g. -c all,0.99 or ITS,0.97,LSU,0.99." << endl;
    #ifdef DATABASE
    std::cout << "-D / --db_type [sqlite/fasta] Sets if the database is in sqlite or fasta format." << endl;
    #endif // DATABASE
    std::cout << "--help / -h                   Print this this help text." << endl;
    std::cout << "--min_length / -m [1+]        Sets the minimum length of sequences to consider for clustering, e.g. -m 100." << endl;
    std::cout << "--no_cluster / -n             Turn clustering off. Only calculating alignment groups." << endl;
    std::cout << "--previous_clusters / -p      Only sequences with cluster marked as 'lead' will be considered for further clustering." << endl;
    #ifdef PTHREAD
    std::cout << "--threads / -T [1+]           Set the number of threads additional to the controling thread, e.g. -T 4." << endl;
    #endif /* PTHREAD */
}

/*** Function to cluster each table ***/
void cluster_each_table ( const char* file, const char* databasetype, const string cut_off, const int min_length, bool only_lead, bool perform_clustering ) {
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
		std::cout << "Could not find appropriate cut off (" << present_cut_off<< ") for " << *table << ". Will only define alignment groups and not cluster." << endl;
		continue;
	    }
	}
	std::cout << "Checking " << *table << endl;
	cluster(database, *table, present_cut_off, min_length, only_lead);
	if (!strcmp(databasetype,"fasta")) {
	    database.print_clusters();
	}
    }
}

void cluster( seqdatabase& db, const string table, const float cut_off, const int min_length, bool only_lead ) {
    #ifdef PTHREAD
    pthread_t thread[n_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    void *status;
    string sequence1;
    string accno1;
    string cluster1;
    sequence_package two_sequences[n_threads];
    int next_thread=0;
    bool activated[n_threads];
    for (int i=0; i < n_threads; ++i) activated[i] = false;
    #else /* PTHREAD */
    sequence_package two_sequences;
    #endif /* PTHREAD */
    align_group deviations;
    char mode='1';
    stringstream converter;
    converter << min_length;
    string char_min_length(converter.str());

    std::cout << "Starting pairwise alignment." << endl;
    string previous = "empty";
    //Need a fasta parser here as alternative, we should prepars
    //Maybe make an object that does the different steps and keep track
    if (db.initiate_sequence_retrieval(table, char_min_length)) {
    	while(!db.all_pairs()) {
    	    db.move_to_next_pair(only_lead);
    	    #ifdef PTHREAD
    	    if (next_thread >= n_threads) next_thread = 0;
    	    int thread_code=0;
    	    if (activated[next_thread]) thread_code = pthread_join(thread[next_thread], &status);
    	    if (thread_code) {
    		std::cout << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
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
    		thread_code = pthread_create(&thread[next_thread], &attr, pthread_align_pair, (void *) &two_sequences[next_thread]);
    		activated[next_thread] = true;
    		++next_thread;
    		if (thread_code) {
    		    std::cout << "ERROR!!! Return code from pthread_create() is: " << thread_code << endl;
    		    exit(-1);
    		}
    	    }
    	    #else /* PTHREAD */
    	    two_sequences.accno1 = db.get_accno1();
    	    two_sequences.sequence1 = db.get_sequence1();
    	    two_sequences.accno2 = db.get_accno2();
    	    two_sequences.sequence2 = db.get_sequence2();
    	    two_sequences.table = &table;
    	    two_sequences.db = &db;
    	    two_sequences.cut_off = &cut_off;
    	    two_sequences.deviations = &deviations;
    	    align_pair ( &two_sequences ); // this is where the action is?
    	    #endif /* PTHREAD */
	    cout << '.';
	}
    }
    else cerr << "Could not initiate sequence retrieval. No aligning done for " << table << "." << endl;
    #ifdef PTHREAD
    pthread_attr_destroy(&attr);
    for (int i=0; i < n_threads; ++i) {
        int thread_code = pthread_join(thread[i], &status);
        if (thread_code) {
            std::cout << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
            exit (-1);
        }
    }
    #endif /* PTHREAD */
    cout << endl;
    std::cout << "Finished aligning. Calculating mad to determine taxonomic level suitable for alignment." << endl;
    cout << "Alignment groups for " << table << endl;
    string alignment_groups = deviations.get_levels( table );
    std::cout << "    Alignment groups: " << alignment_groups << endl;
    #if DATABASE
    db.insert_alignment_group(table, alignment_groups); 
    #endif //DATABASE
    #ifdef DEBUG
    cout << "    Aprox. mad. entire group = " << deviations.aprox_mad() << endl;
    deviations.print_hierarchy();
    #endif /* DEBUG */
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
    sequences.set_cost_matrix( 7, -5 );
    sequences.align();
    #ifdef PTHREAD
    pthread_mutex_lock (&databasemutex);
    #endif /* PTHREAD */
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
    #ifdef PTHREAD
    pthread_mutex_unlock(&databasemutex);
    #endif /* PTHREAD */
}




