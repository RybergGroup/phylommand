#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "seqpair.h"
#include "sqlite3.h"
#include "align_group.h"
#ifdef PTHREAD // if compiling multithreads
#include <pthread.h>
#endif /*PTHREAD*/

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
    //float proportion_N1;
    //float proportion_N2;
    const string *table;
    sqlite3 *db;
    const float *cut_off;
    align_group *deviations;
};

#ifdef PTHREAD
void *pthread_align_pair (void *threadarg);
#endif /* PTHREAD */

void align_pair ( sequence_package *two_sequences );
void clust_update( string accno, string cluster, string table, sqlite3 *db, bool where_accno);
void cluster( sqlite3 *db, const string table, const float cut_off, const int min_length, bool only_lead );
string get_cluster( string accno, string table, sqlite3 *db);
float get_comp_value ( const string accno, const string table, sqlite3 *db);
void cluster_each_table ( const char* database, const string cut_off, const int min_length, bool only_lead, bool perform_clustering );
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
   // If cut off value out of bound, give instructions then  quit
/*   if (cut_off <= 0.0 || cut_off > 1.1) {
       std::cout << "Cut off (--cut_off or -c) must be a value between 0 and 1.1, e.g. 0.99." << endl;
       error_flag = 1;
   }*/
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
   // If database can not be opened, give instructions then quit
   std::cout << "Using database: " << argv[argc-1] << '.' << endl;
   sqlite3 *db;
   if ( sqlite3_open(argv[argc-1], &db) != 0 ) {
       std::cout << "Could not open " << argv[argc-1] << ". Last argument has to be a sqlite database." << endl;
       error_flag = 1;
   }
   sqlite3_close(db);
   // if any problems above, quit!
   if (error_flag) {
       std::cout << "Quitting quietly." << endl;
       return 0;
   }

   std::cout << "Started at:" << endl << asctime( timeinfo );
   // Execute clustering, last argument is interpreted as the database file
   cluster_each_table ( argv[argc-1], cut_off, min_length, only_lead, perform_clustering );
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
    std::cout << "--cut_off / -c [0-1]        Sets the cut off value in pairwise similarity for clustering of each gene. Should be" << endl;
    std::cout << "                                comma separated string with gene name first and value second. If all genes should" << endl;
    std::cout << "                                have same cut off 'all' could be given as gene. E.g. -c all,0.99 or ITS,0.97,LSU,0.99." << endl;
    std::cout << "--help / -h                 Print this this help text." << endl;
    std::cout << "--min_length / -m [1+]      Sets the minimum length of sequences to consider for clustering, e.g. -m 100." << endl;
    std::cout << "--no_cluster / -n           Turn clustering off. Only calculating alignment groups." << endl;
    std::cout << "--previous_clusters / -p    Only sequences with cluster marked as 'lead' will be considered for further clustering." << endl;
    #ifdef PTHREAD
    std::cout << "--threads / -T [1+]         Set the number of threads additional to the controling thread, e.g. -T 4." << endl;
    #endif /* PTHREAD */
}

/*** Function to cluster each table ***/
void cluster_each_table ( const char* database, const string cut_off, const int min_length, bool only_lead, bool perform_clustering ) {
    sqlite3 *db;
    if ( sqlite3_open(database, &db) == 0 ) {
        sqlite3_stmt *statement;
        string query = "CREATE TABLE alignment_groups (gene TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', tree TEXT DEFAULT 'empty', tree_method TEXT DEFAULT 'empty', alignable INTEGER, PRIMARY KEY (gene, taxon));";
        if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
            if (sqlite3_step(statement) != SQLITE_DONE) {
                std::cerr << "Could not create table for alignment_groups (1). Quitting." << endl;
                return;
            }
            sqlite3_finalize(statement);
        }
        else {
            sqlite3_finalize(statement);
            query = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;";
            bool flag = 0;
            if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
                while (sqlite3_step(statement) == SQLITE_ROW) {
                    string table = (char*)sqlite3_column_text(statement,0);
                    if (!table.compare("alignment_groups")) {
                        flag = 1;
                        break;
                    }
                }
            }
            sqlite3_finalize(statement);
            if (!flag) {
                std::cerr << "Could not create table for alignment_groups (2). Quitting." << endl;
                return;
            }
        }
        query = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;";
        if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
            while (sqlite3_step(statement) == SQLITE_ROW) {
                string table = (char*)sqlite3_column_text(statement,0);
                if (!table.compare("gb_data") || !table.compare("alignments") || !table.compare("alignment_groups")) continue;
                float present_cut_off=0;
                if (perform_clustering) {
                    int length = cut_off.length();
                    string gene;
                    for (int i=0; i < length; ++i) {
                        if (cut_off[i]==',') {
                            if(!table.compare(gene) || !table.compare("all")) {
                                string number;
                                ++i;
                                while (cut_off[i]!=',' && i<length) {
                                    number+=cut_off[i];
                                    ++i;
                                }
                                present_cut_off = atof(number.c_str());
                                break;
                            }
                        }
                        else gene+=cut_off[i];
                    }
                    if (present_cut_off < 0.000000001) {
                        std::cout << "Could not find appropriate cut off for " << table << ". Will only define alignment groups and not cluster." << endl;
                        continue;
                    }
                }
                std::cout << "Checking " << table << endl;
                cluster(db, table, present_cut_off, min_length, only_lead);
            }
        }
        sqlite3_finalize(statement);
    }
    sqlite3_close(db);
}

void cluster( sqlite3 *db, const string table, const float cut_off, const int min_length, bool only_lead ) {
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
    char char_min_length[33];
    sprintf (char_min_length, "%d", min_length);

    std::cout << "Starting pairwise alignment." << endl;
    string query = "SELECT accno,sequence,cluster FROM ";
    query += table;
    query += " WHERE LENGTH(sequence)>";
    query += char_min_length;
    query += ';';
    string previous = "empty";
    while(mode!='9') {
        sqlite3_stmt *statement;
        if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
            if (previous.compare("empty")) { // if previous is not empty
                while (1) {
                    string accno;
                    if (sqlite3_step(statement) == SQLITE_ROW) accno = (char*)sqlite3_column_text(statement,0);
                    else {
                        mode = '9';
                        break;
                    }
                    if (!accno.compare(previous)) break;
                }
            }
            while (1) { 
                if (sqlite3_step(statement) == SQLITE_ROW) {
                    if (only_lead) {
                        string cluster = (char*)sqlite3_column_text(statement,2);
                        if (cluster.compare("lead")) {
                            if (mode == '1') previous = (char*)sqlite3_column_text(statement,0); // if looking for first seq save accno
                            continue; // if not lead sequence continue
                        }
                    }
                    if (mode == '1') {
                        #ifdef PTHREAD
                        accno1 = (char*)sqlite3_column_text(statement,0);
                        previous = accno1;
                        sequence1 = (char*)sqlite3_column_text(statement,1);
                        #else /* PTHREAD */
                        two_sequences.accno1 = (char*)sqlite3_column_text(statement,0);
                        two_sequences.sequence1 = (char*)sqlite3_column_text(statement,1);
                         previous = two_sequences.accno1;
                        #endif /* PTHREAD */
                        mode = '2';
                    }
                    else if ( mode == '2' || mode == '3') {
                        #ifdef PTHREAD
                        if (next_thread >= n_threads) next_thread = 0;
                        int thread_code=0;
                        if (activated[next_thread]) thread_code = pthread_join(thread[next_thread], &status);
                        if (thread_code) {
                            std::cout << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
                            exit (-1);
                        }
                        else {
                            two_sequences[next_thread].accno1 = accno1;
                            two_sequences[next_thread].sequence1 = sequence1;
                            two_sequences[next_thread].accno2 = (char*)sqlite3_column_text(statement,0);
                            two_sequences[next_thread].sequence2 = (char*)sqlite3_column_text(statement,1);
                            two_sequences[next_thread].table = &table;
                            two_sequences[next_thread].db = db;
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
                        two_sequences.accno2 = (char*)sqlite3_column_text(statement,0);
                        two_sequences.sequence2 = (char*)sqlite3_column_text(statement,1);
                        two_sequences.table = &table;
                        two_sequences.db = db;
                        two_sequences.cut_off = &cut_off;
                        two_sequences.deviations = &deviations;
                        align_pair ( &two_sequences );
                        #endif /* PTHREAD */
                        if (mode == '2') {
                            mode = '3';
                        }
                    }
                    #ifdef PTHREAD
                    #else /* PTHREAD */
                    two_sequences.accno2.clear();
                    two_sequences.sequence2.clear();
                    #endif /* PTHREAD */
                }
                else {
                    if ( mode == '2' || mode == '1') mode = '9';
                    break;
                }
            }
            #ifdef PTHREAD
            #else /* PTHREAD */
            two_sequences.accno1.clear();
            two_sequences.sequence1.clear();
            #endif /* PTHREAD */
            if (mode !='9') mode='1';
        }
        else {
            std::cerr << "Could not prepare SQLite query!!! Quitting." << endl;
            exit(-1);
        }
        sqlite3_finalize(statement);
    }
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
    std::cout << "Finished aligning. Calculating mad to determine taxonomic level suitable for alignment." << endl;
    cout << "Alignment groups for " << table << endl;
    string alignment_groups = deviations.get_levels( db, table );
    std::cout << "    Alignment groups: " << alignment_groups << endl;
    if (alignment_groups.compare("empty")) { // if something in alignment groups
        string taxon;
        char alignable;
        for (int i=0; i < alignment_groups.length(); ++i) {
            if (alignment_groups[i] == ';') {
                sqlite3_stmt *statement;
                string insert = "INSERT INTO alignment_groups (gene, taxon, alignable) VALUES ('";
                insert += table;
                insert += "', '";
                insert += taxon;
                insert += "', ";
                insert += alignable;
                insert += ");";
                if(sqlite3_prepare_v2(db, insert.c_str(), -1, &statement, 0) == SQLITE_OK) {
                    sqlite3_step(statement);
                }
                else std::cerr << "Failed to uppdate database. SQL statement: " << insert << endl << "Proceeding reluctantly." << endl;
                sqlite3_finalize(statement);
                std::cout << "    " << taxon << "    ";
                if (alignable == '1') std::cout << "Alignable" << endl;
                else std::cout << "Not alignable" << endl;
                taxon.clear();
            }
            else if (alignment_groups[i] == '_') {
                if (alignment_groups[++i] == 'A') alignable = '1';
                else alignable = '0';
            }
            else taxon += alignment_groups[i];
        }
    }
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
    sequences.set_cost_matrix( 7, -5 );
    sequences.align();
    #ifdef PTHREAD
    pthread_mutex_lock (&databasemutex);
    #endif /* PTHREAD */
    if (*two_sequences->cut_off > 0.000000001) {
        string cluster1 = get_cluster( two_sequences->accno1, *two_sequences->table, two_sequences->db );
        string cluster2 = get_cluster( two_sequences->accno2, *two_sequences->table, two_sequences->db );
        if (sequences.similarity(1) > *two_sequences->cut_off) {
            float comp1;
            float comp2;
            if (!cluster1.compare( "empty" ) || !cluster1.compare( "lead" )) comp1 = get_comp_value( two_sequences->accno1, *two_sequences->table, two_sequences->db );
            else comp1 = get_comp_value( cluster1, *two_sequences->table, two_sequences->db );
            if (!cluster2.compare( "empty" ) || !cluster2.compare( "lead" )) comp2 = get_comp_value( two_sequences->accno2, *two_sequences->table, two_sequences->db );
            else comp2 = get_comp_value( cluster2, *two_sequences->table, two_sequences->db );
            // If sequence 1 better
            if (comp1 >= comp2) {
                if ( !cluster1.compare( "empty" ) ) {
                    clust_update( two_sequences->accno1, "lead", *two_sequences->table, two_sequences->db, 1 );
                    clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 1 );
                    if (!cluster2.compare( "lead" ) || !cluster2.compare( "empty" )) clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 0 );
                    else {
                        clust_update( cluster2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 1 );
                        clust_update( cluster2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 0 );
                    }
                }
                else if ( !cluster1.compare( "lead" ) ) {
                    clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 1 );
                    if (!cluster2.compare( "lead" ) || !cluster2.compare( "empty" )) clust_update( two_sequences->accno2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 0 );
                    else {
                        clust_update( cluster2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 1 );
                        clust_update( cluster2, two_sequences->accno1, *two_sequences->table, two_sequences->db, 0 );
                    }
                }
                else {
                    clust_update( two_sequences->accno2, cluster1, *two_sequences->table, two_sequences->db, 1 );
                    if (!cluster2.compare( "lead" ) || !cluster2.compare( "empty" )) clust_update( two_sequences->accno2, cluster1, *two_sequences->table, two_sequences->db, 0 );
                    else {
                        clust_update( cluster2, cluster1, *two_sequences->table, two_sequences->db, 1 );
                        clust_update( cluster2, cluster1, *two_sequences->table, two_sequences->db, 0 );
                    }
                }
            }
         // If sequence 2 better
            else {
                if ( !cluster2.compare( "empty" ) ) { // if no previous annotation
                    clust_update( two_sequences->accno2, "lead", *two_sequences->table, two_sequences->db, 1 ); // set best sequence to lead
                    clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 1 ); // set worse sequence to the same cluster
                    if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 0 ); // if worse sequence lead or empty uppdate sequences in that cluster
                    else { // if not lead
                        clust_update( cluster1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 1 ); // change lead of cluster to point to new best sequence
                        clust_update( cluster1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 0 ); // change all sequences previusly pointing to that sequence to point to new best sequence
                    }
                }
                else if ( !cluster2.compare( "lead" ) ) {
                    clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 1 ); // set worse sequence to the same cluster
                    if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) clust_update( two_sequences->accno1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 0 ); // if worse sequence lead or empty uppdate sequences in that cluster
                    else {
                        clust_update( cluster1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 1 ); // change lead of cluster to point to new best sequence
                        clust_update( cluster1, two_sequences->accno2, *two_sequences->table, two_sequences->db, 0 ); // change all sequences previusly pointing to that sequence to point to new best sequence
                    }
                }
                else {
                    clust_update( two_sequences->accno1, cluster2, *two_sequences->table, two_sequences->db, 1 ); // set worse sequence to point to the same sequence as better sequence
                    if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) clust_update( two_sequences->accno1, cluster2, *two_sequences->table, two_sequences->db, 0 ); // set worse sequence to point to same sequence as better sequence
                    else {
                        clust_update( cluster1, cluster2, *two_sequences->table, two_sequences->db, 1 ); // change lead of worse cluster to point to the same sequence as new best sequence
                        clust_update( cluster1, cluster2, *two_sequences->table, two_sequences->db, 0 ); // change rest of cluster to point to the same sequence as new best sequence
                    }
                }
            }
        }
        else {
            if ( !cluster1.compare( "empty" ) ) clust_update( two_sequences->accno1, "lead", *two_sequences->table, two_sequences->db, 1 );
            if ( !cluster2.compare( "empty" ) ) clust_update( two_sequences->accno2, "lead", *two_sequences->table, two_sequences->db, 1 );
            if ((!cluster1.compare( "empty" ) || !cluster1.compare( "lead" )) && (!cluster2.compare( "empty" ) || !cluster2.compare( "lead" ))) {
                two_sequences->deviations->insert_value( two_sequences->accno1, two_sequences->accno2, sequences.jc_distance()-(1-sequences.similarity()), two_sequences->db );
            }
        }
    }
    else two_sequences->deviations->insert_value( two_sequences->accno1, two_sequences->accno2, sequences.jc_distance()-(1-sequences.similarity()), two_sequences->db );
    #ifdef PTHREAD
    pthread_mutex_unlock(&databasemutex);
    #endif /* PTHREAD */
}

void clust_update( const string accno, const string cluster, const string table, sqlite3 *db, const bool where_accno) {
    sqlite3_stmt *statement;
    string update = "UPDATE ";
    update += table;
    update += " SET cluster='";
    update += cluster;
    if (where_accno) {
        update += "' WHERE accno='";
        update += accno;
    }
    else {
        update += "' WHERE cluster='";
        update += accno;
    }
    update += "';";
    if(sqlite3_prepare_v2(db, update.c_str(), -1, &statement, 0) == SQLITE_OK) {
        sqlite3_step(statement);
    }
    sqlite3_finalize(statement);
}

string get_cluster( const string accno, const string table, sqlite3 *db) {
    string cluster;
    sqlite3_stmt *statement;
    string query = "SELECT cluster FROM ";
    query += table;
    query += " WHERE accno='";
    query += accno;
    query += "';";
    if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
        if (sqlite3_step(statement) == SQLITE_ROW) {
            cluster=(char*)sqlite3_column_text(statement,0);
        }
    }
    sqlite3_finalize(statement);
    return cluster;
}

float get_comp_value ( const string accno, const string table, sqlite3 *db) {
    int length=0;
    float prop_N=0.0;
    string query  = "SELECT LENGTH(";
    query += table;
    query += ".sequence),gb_data.proportion_N FROM gb_data INNER JOIN ";
    query += table;
    query += " on gb_data.accno=";
    query += table;
    query += ".accno WHERE gb_data.accno='";
    query += accno;
    query += "';";
    sqlite3_stmt *statement;
    if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
        if (sqlite3_step(statement) == SQLITE_ROW) {
            length=sqlite3_column_int(statement,0);
            prop_N=sqlite3_column_double(statement,1);
        }
    }
    return length*(1-prop_N);
}
