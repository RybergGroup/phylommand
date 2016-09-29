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
#include "argv_parser.h"

using namespace std;

unsigned int n_threads=1;
#ifdef PTHREAD
pthread_mutex_t databasemutex;
class print_queue {
    public:
	print_queue():start(0) {};
	~print_queue() {
	    delete_queue(start);
	}
	stringstream* new_position ();
	string print_next ( ) { return print_next( false ); };
	string print_next( bool print_empty );
	bool empty() {
	    if (start == 0) return true;
	    else return false;
	};
	bool something_in_first_position() {
	    if (start!=0 && start->output.rdbuf()->in_avail()) return true;
	    else return false;
	};
	void pop();
    private:
    class node {
	public:
	node (): next(0) {};
	stringstream output;
	node* next;
    };
    node* start;
    void delete_queue ( node* position ); 
    node* go_to_end( node* position ) {
	if (position == 0) return 0;
	else if (position->next != 0) return go_to_end(position->next);
	else return position;
    };
};
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
    ostream* outputstream;
};

#ifdef PTHREAD
void *pthread_align_pair (void *threadarg);
#endif /* PTHREAD */

void align_pair ( sequence_package *two_sequences );
void cluster ( seqdatabase& db, const string table, const float cut_off, const int min_length, const bool only_lead, const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet );
void cluster_each_table ( const string& file, const char* databasetype, const string cut_off, const int min_length, const bool only_lead, /*const bool perform_clustering,*/ const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet, const string& taxonomy_file );

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
    string taxonomy_file_name;
    ifstream infile;
    stringstream stdin_holder;
    // From alignment_groups
    time_t rawtime = time(0);
    struct tm * timeinfo = localtime( &rawtime );
    // Cut off value set to default
    string cut_off("all,0.999");
    int min_length(100);
    bool only_lead(false);
    //bool perform_clustering = true;
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
	else if (!strcmp(argv[i], "--format")) {
	    if (i+1 < argc && argv[i+1][0] != '-') {
		++i;
		if (!strcmp(argv[i],"fasta")) strcpy(db_type, "fasta");
		else if (!strcmp(argv[i],"pairfst") || !strcmp(argv[i],"pairfa")) strcpy(db_type, "pairfa");
		#ifdef DATABASE
		else if (!strcmp(argv[i],"sqlite")) strcpy(db_type, "sqlite");
		#endif //DATABASE
		else {
		    cerr << argv[i] << " is not a valid option for --format. The options are ";
		    #ifdef DATABASE
		    cerr << "sqlite, ";
		    #endif //DATABASE
		    cerr << "pairfst, or fasta." << endl;
		    return 1;
		}
	    }
	    else {
		cerr << "--format require fasta";
		#ifdef DATABASE
		cerr << ", sqlite,";
		#endif //DATABASE
		cerr << " or pairwise as extra argument. Use -h for more help." << endl;
		return 1;
	    }
	}
	else if ( !strcmp(argv[i],"-g") || !strcmp(argv[i],"--group")) {
	    output_mode = 'A';
	    if (i+1 < argc && argv[i+1][0] != '-') {
		++i;
		vector<string> arguments;
		argv_parser::pars_sub_args(argv[i], ':', arguments);
		#ifdef DEBUG
		cerr << "First extra argument: " << arguments[0] << endl;
		#endif //DEBUG
		if (!arguments[0].compare("alignment_groups")) output_mode = 'A';
		else if (!arguments[0].compare("cluster")) output_mode = 'C';
		else if (!arguments[0].compare("both")) output_mode = 'B';
		else {
		    cerr << "Do not recognize argument '" << arguments[0] << "'. Options are 'alignment_groups', 'cluster', 'both'." << endl;
		    return 1;
		}
		if (arguments.size() > 1) {
		    for (vector<string>::iterator i=arguments.begin()+1; i != arguments.end(); ++i) {
			vector<string> sub_args;
		       	argv_parser::pars_sub_args(i->c_str(),'=', sub_args);
			if ((!sub_args[0].compare("cut-off") || !sub_args[0].compare("cut_off") || !sub_args[0].compare("cutoff")
			    || !sub_args[0].compare("cut off")) && sub_args[0].size() > 1) cut_off = sub_args[1];
			else if ((!sub_args[0].compare("min-length") || !sub_args[0].compare("min_length") || !sub_args[0].compare("minlength")
			    || !sub_args[0].compare("min length")) && sub_args[0].size() > 1) min_length = atoi(sub_args[1].c_str());
			else if ((!sub_args[0].compare("taxon") || !sub_args[0].compare("taxonomy") || !sub_args[0].compare("taxon_file")
			    || !sub_args[0].compare("taxonfile")) && sub_args[0].size() > 1) taxonomy_file_name = sub_args[1];
			else if (!sub_args[0].compare("only_lead") || !sub_args[0].compare("only-lead") || !sub_args[0].compare("only lead")
			    || !sub_args[0].compare("previous_clusters")) only_lead = 1;
			else {
			    cerr << "Unrecognized argument or missing value for argument given to -g/--group: " << sub_args[0] << "." << endl;
			    return 1;
			}
		    }
		}
		#ifdef DEBUG
		cerr << "Output: " << output_mode << " Cut-off: " << cut_off;
		#ifdef DATABASE
		cerr << " Database type: " << db_type;
		#endif //DATABASE
		cerr << endl;
		#endif //DEBUG
	    }
	}
	/// From alignment_groups
	/*else if ( !strcmp(argv[i],"-c") || !strcmp(argv[i],"--cut_off") ) {
	    ++i;
	    if (i < argc && argv[i][0] != '-') { 
		cut_off = argv[i];
	    }
	    else {
	       	std::cerr << "--cut_off or -c must be followed by a number or a comma separated string. Quitting quietly." << endl;
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
	}*/
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
	cluster_each_table ( file_name, db_type, cut_off, min_length, only_lead, /*perform_clustering,*/ aligned, output_mode, output_names, matrix, quiet, taxonomy_file_name );
	// Calculate end time
	rawtime = time(0);
	timeinfo = localtime( &rawtime );
	if (!quiet) std::cerr << "Ended at:" << endl << asctime( timeinfo );
     	return 0; // return normally
}

/*** Function to print help ***/
void help() {
    std::cout << "Pairalign " << VERSION << " will perform pairwise alignment of DNA sequences given in fasta" << endl; 
    std::cout << "format through standard in." << endl;
    std::cout << "(c) Martin Ryberg " << YEAR << "." << endl << endl;
    std::cout << "Usage:" << endl << "pairalign [arguments] < inputfile.fasta" << endl << "pairalign [arguments] inputfile.fasta" << endl << endl;
    std::cout << "Arguments:" << endl;
    std::cout << "--aligned / -A                  input file is already aligned." << endl;
    std::cout << "--alignments / -a               output aligned sequences pairwise." << endl;
////////////////////////////
    /*std::cout << "--cut_off / -c [value/s]        sets the cut off value in pairwise similarity" << endl;
    std::cout << "                                for clustering of each gene. Should be comma" << endl;
    std::cout << "                                separated string with gene name first and value" << endl;
    std::cout << "                                second. If all genes should have same cut off" << endl;
    std::cout << "                                'all' could be given as gene. E.g. -c all,0.99" << endl;
    std::cout << "                                or ITS,0.97,LSU,0.99." << endl;
    #ifdef DATABASE
    std::cout << "--db_type / -D [sqlite/fasta]   sets if the database is in sqlite or fasta" << endl;
    std::cout << "                                format." << endl;
    #endif // DATABASE*/
/////////////////////////////////
    std::cout << "--difference / -i               output difference between the Jukes-Cantor (JC)" << endl;
    std::cout << "                                distance and proportion different sites." << endl;
    std::cout << "--distances / -d                output proportion different sites, JC distance," << endl;
    std::cout << "                                and diference between the two." << endl;
    #ifdef DATABASE
    std::cout << "--format [fasta/pairfst/sqlite] ";
    #else
    std::cout << "--format [fasta/pairfst]        ";
    #endif //DATABASE
    std::cout << "set the format of the input to fasta or fasta" << endl; // add that fasta is default
    std::cout << "                                with sequences pairwise (as output given the -a" << endl; 
    std::cout << "                                -n option). If sequences are aligned give the -A" << endl;
    std::cout << "                                switch.";
    #ifdef DATABASE
    std::cout << " If the sqlite option is given, input" << endl;
    std::cout << "                                is expected to be a SQLite database file with" << endl;
    std::cout << "                                one table named gb_data that have a column named" << endl;
    std::cout << "                                accno containing accession number for each" << endl; 
    std::cout << "                                sequence, and one column named taxon_string" << endl;
    std::cout << "                                containing the taxonomic hierarchy for the" << endl;
    std::cout << "                                sequence, in the format given for --group. It" << endl;
    std::cout << "                                should also have a table for each gene region to" << endl;
    std::cout << "                                be tested. Each of these tables should have a" << endl;
    std::cout << "                                column named accno with the accession for each" << endl;
    std::cout << "                                sequence, one column named sequence with the" << endl;
    std::cout << "                                sequence, and one column named cluster (with" << endl;
    std::cout << "                                empty as default value).";
    #endif //DATABASE
    std::cout << endl;
    std::cout << "--group / -g                    this option will cluster sequences that are" << endl;
    std::cout << "                                similar and/or find the most inclusive taxa in a" << endl; 
    std::cout << "                                hierarchy that are alignable according to MAD" << endl;
    std::cout << "                                (Smith et al. 2009, BMC evol. Biol. 9:37). It" << endl;
    std::cout << "                                need the taxonomy given after a (the first) | in" << endl;
    std::cout << "                                the sequence name or in a separate file. Each" << endl;
    std::cout << "                                taxa in the hierarchy should be separated by a" << endl;
    std::cout << "                                semicolon, with the highest rank first and then" << endl;
    std::cout << "                                increasingly nested levels until the lowest" << endl;
    std::cout << "                                known level for the sequence. The groups that" << endl;
    std::cout << "                                can be aligned are put in a file with the ending" << endl;
    std::cout << "                                .alignment_groups and printed to the screen" << endl;
    std::cout << "                                preceded by #. Clusters are printed to the" << endl;
    std::cout << "                                screen after a heading, preceded by ###. To get" << endl; 
    std::cout << "                                alignable groups give 'alignment_groups' as" << endl;
    std::cout << "                                extra argument, to cluster give 'cluster', and" << endl;
    std::cout << "                                to do both give 'both'. Cut off value for" << endl;
    std::cout << "                                pairwise similarity can be given after colon (:)" << endl;
    std::cout << "                                by cut-off= followed value, e.g. -g both:" << endl;
    std::cout << "                                cut-off=0.97. A file with taxonomy can be given" << endl;
    std::cout << "                                with taxonomy=. The taxonomy file should have" << endl;
    std::cout << "                                the taxonomy (as above) first on each row" << endl;
    std::cout << "                                followed by a |, and the sequence name with that" << endl;
    std::cout << "                                taxonomy as a comma (,) and/or space ( )" << endl;
    std::cout << "                                separated string. The same taxon can be repeated" << endl;
    std::cout << "                                several times.";
    #ifdef DATABASE
    std::cout << " When using the sqlite format" << endl;
    std::cout << "                                'only_lead' can be added as an extra argument" << endl;
    std::cout << "                                behind a colon to only cluster sequences marked" << endl;
    std::cout << "                                as 'lead' in the cluster column and min_length=" << endl;
    std::cout << "                                followed by length of sequences to only consider" << endl;
    std::cout << "                                sequences longer than given length. E.g. -g both" << endl;
    std::cout << "                                :only_lead:db_type=sqlite:cut-off=0.97:" << endl;
    std::cout << "                                min_length=100. Output will be saved to the" << endl;
    std::cout << "                                database.";
    #endif //DATABASE
    std::cout << endl;
    std::cout << "--help / -h                     print this help." << endl;
    std::cout << "--jc_distance / -j              output Jukes-Cantor (JC) distance." << endl;
    std::cout << "--matrix / -m                   output in the form of a space separated" << endl;
    std::cout << "                                left-upper triangular matrix." << endl;
    //std::cout << "--min_length / -m [value]       sets the minimum length of sequences to consider"<< endl;
    //std::cout << "                                for clustering, e.g. -m 100." << endl;
    std::cout << "--names / -n                    output sequence names (if outputting alignments" << endl;
    std::cout << "                                then in fasta format)." << endl;
    //std::cout << "--no_cluster / -O               turn clustering off when running" << endl;
    //std::cout << "                                --alignment_groups. Only calculating alignment" << endl;
    //std::cout << "                                groups." << endl;
    //#ifdef DATABASE
    //std::cout << "--previous_clusters / -P        only consider sequences with cluster marked as" << endl;
    //std::cout << "                                'lead' for further clustering." << endl;
    //#endif // DATABASE
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
void cluster_each_table ( const string& file, const char* databasetype, const string cut_off, const int min_length, const bool only_lead, /*const bool perform_clustering,*/ const bool aligned, const char output, const bool output_names, const bool matrix, const bool quiet, const string& taxonomy_file ) {
    seqdatabase database(file,databasetype);
    if (!quiet) cerr << "Opened " << databasetype << " database." << endl;
    if (!taxonomy_file.empty()) {
	if (!quiet) cerr << "Parsing taxonomy from " << taxonomy_file << "." << endl;
	ifstream file(taxonomy_file.c_str());
	if (file.good()) {
	    map<string,string> taxonomy_map;
	    string taxa;
	    string accno;
	    bool taxonomy(true);
	    while(file) {
		char input_char = file.get();
		if (input_char == '|' && taxonomy) taxonomy = false;
		else if (!taxonomy && (input_char == ' ' || input_char == ',' || input_char == '\n' || input_char == '\r') && !accno.empty() && !taxa.empty()) {
		    taxonomy_map[accno] = taxa;
		    #ifdef DEBUG
		    cerr << "Added " << taxa << " for " << accno << "." << endl;
		    #endif //DEBUG
		    accno.clear();
		    if (input_char == '\r' || input_char == '\n') { taxonomy = true; taxa.clear(); }
		}
		else if (input_char == '\r' || input_char == '\n') { taxonomy = true; taxa.clear(); }
		else if (taxonomy) taxa += input_char;
	       	else if (!taxonomy) accno += input_char;
		
	    }
	    if (!taxonomy_map.empty()) {
		database.add_taxonomy(taxonomy_map);
		if (!quiet) cerr << "Added taxonomy." << endl;
	    }
	    else cerr << "Was not able to pars taxonomy from " << taxonomy_file << "." << endl;
	}
	else cerr << "Was not able to open " << taxonomy_file << ". No taxonomy read." << endl;
    }
    else if (!strcmp(databasetype,"pairfa")) {
	map<string,string> taxonomy_map;
	taxonomy_map["default"] = "all";
	database.add_taxonomy(taxonomy_map);
	if (!quiet && (output == 'A' || output == 'B')) cerr << "All sequences will be treated as from same taxon." << endl;
    }
    if ((output == 'A' || output == 'B') && !database.alignment_groups_present()) {
	if (!quiet) cerr << "No alignment_groups file/table present. Trying to create it." << endl;
	if (!database.create_alignment_groups()) {
	    cerr << "Failed to create alignment_groups." << endl;
	    return;
	}
    }
    vector<string> tables = database.tables_in_database();
    for (vector<string>::const_iterator table = tables.begin(); table != tables.end(); ++table) {
	if (!table->compare("gb_data") || !table->compare("alignments") || !table->compare("alignment_groups")) continue;
	float present_cut_off=0.0;
	if (output == 'B' || output == 'C') { //perform_clustering) {
	    int length = cut_off.length();
	    string gene;
	    int i(0);
	    for (; i < length; ++i) {
		if (cut_off[i]==',') {
		    if(!gene.compare(*table) || !gene.compare("all")) {
			string number;
			++i;
			while (i<length && cut_off[i]!=',') {
			    number += cut_off[i];
			    ++i;
			}
			present_cut_off = atof(number.c_str());
			gene.clear();
			break;
		    }
		}
		else gene+=cut_off[i];
	    }
	    if (i > 0 && i >= length && !gene.empty() && (gene[0] == '0' || gene[0] == '1' || gene[0] == '2' || gene[0] == '3' ||
		    gene[0] == '4' || gene[0] == '5' || gene[0] == '6' || gene[0] == '7' || gene[0] == '8' || gene[0] == '9' ||
		    gene[0] == '.') ) present_cut_off = atof(gene.c_str());
	    if (present_cut_off < 0.000000001) {
		std::cerr << "Could not find appropriate cut off (" << present_cut_off<< ") for " << *table << ". Will only define alignment groups and not cluster." << endl;
		continue;
	    }
	    if (!quiet) cerr << "Using the cut off: " << present_cut_off << "." << endl;
	}
	if (!quiet) std::cerr << "Checking " << *table << endl;
	cluster(database, *table, present_cut_off, min_length, only_lead, aligned, output, output_names, matrix, quiet);
	//#ifdef DEBUG
	//cerr << "Output: " << output << " | " << "Database type: " << databasetype << endl;
	//database.print_clusters(cerr);
	//#endif //DEBUG
	if ((output == 'B' || output == 'C') && (!strcmp(databasetype,"fasta") || !strcmp(databasetype,"pairfa"))) {
	    cout << "### Clusters " << *table << ", cut-off: " << present_cut_off << " ###" << endl;
	    database.print_clusters(cout);
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
    unsigned int next_thread=0;
    bool activated[n_threads];
    for (unsigned int i=0; i < n_threads; ++i) activated[i] = false;
    print_queue output_queue;
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
    	    #ifdef PTHREAD
	    if (n_threads > 1) {
		pthread_mutex_lock(&databasemutex);
		if (matrix && new_row) {
		    ++n_seq;
		    stringstream* outputstream = output_queue.new_position();
		    if (n_seq > 1)  *outputstream << endl;
		    if (output_names) *outputstream << db.get_accno1() << ' ';
		    for (unsigned int i=1; i<n_seq; ++i) *outputstream << ' ';
		}
		/*while (output_queue.something_in_first_position()) {
		    cout << output_queue.print_next();
		    #ifdef DEBUG
		    cerr << "Printing" << endl;
		    #endif //DEBUG
		}*/
		pthread_mutex_unlock(&databasemutex);
		if (next_thread >= n_threads) next_thread = 0;
		int thread_code=0;
		if (activated[next_thread]) thread_code = pthread_join(thread[next_thread], &status);
		if (thread_code) {
		    std::cerr << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
		    exit (-1);
		}
		else {
		    #ifdef DEBUG
		    cerr << "Prepairing pthread alignment: ";
		    #endif //DEBUG
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
		    if (output != 'A' && output != 'B' && output != 'C') {
			pthread_mutex_lock(&databasemutex);
			two_sequences[next_thread].outputstream = output_queue.new_position();
			pthread_mutex_unlock(&databasemutex);
		    }
		    else two_sequences[next_thread].outputstream = 0;
		    #ifdef DEBUG
		    cerr << two_sequences[next_thread].accno1 << " " << two_sequences[next_thread].accno2 << endl;
		    #endif //DEBUG
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
		if (matrix && new_row) {
		    ++n_seq;
		    if (n_seq > 1) std::cout << endl;
		    if (output_names) std::cout << db.get_accno1() << ' ';
		    for (unsigned int i=1; i<n_seq; ++i) std::cout << ' ';
		}
		#ifdef DEBUG
		cerr << "Prepairing alignment: ";
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
		two_sequences[0].outputstream = &cout;
		#ifdef DEBUG
		cerr << two_sequences[0].accno1 << " " << two_sequences[0].accno2 << endl;
		#endif //DEBUG
		align_pair ( &two_sequences[0] ); // this is where the action is?
		//#endif /* PTHREAD */
		#ifdef PTHREAD
	    }
    	    #endif /* PTHREAD */
	    if (!quiet) cerr << '.';
	}
    }
    else cerr << "Could not initiate sequence retrieval. No aligning done for " << table << "." << endl;
    #ifdef PTHREAD
    pthread_attr_destroy(&attr);
    if (n_threads > 1) {
	for (unsigned int i=0; i < n_threads; ++i) {
	    int thread_code = pthread_join(thread[i], &status);
	    if (thread_code) {
		std::cerr << "ERROR!!! Return code from pthread_join() is: " << thread_code << endl;
		exit (-1);
	    }
	}
    }
    #ifdef PTHREAD
    while (!output_queue.empty()) cout << output_queue.print_next(true);
    #endif /* PTHREAD */
    if ( matrix && output_names ) std::cout << endl << db.get_accno2() << endl; // from pairalign
    #endif /* PTHREAD */
    if (!quiet) cerr << endl;
    if (output == 'A' || output == 'B') {
	if (!quiet) std::cerr << "Finished aligning. Calculating mad to determine taxonomic level suitable for alignment." << endl;
	cout << "# Alignment groups for " << table << endl;
	string alignment_groups = deviations.get_levels( table );
	std::cout << "#    Alignment groups: " << alignment_groups << endl;
	if (!quiet) cerr << "    Aprox. mad. entire group = " << deviations.aprox_mad() << endl;
	//#if DATABASE
	db.insert_alignment_group(table, alignment_groups); 
	//#endif //DATABASE
	#ifdef DEBUG
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
    cerr << "Starting to process: " << two_sequences->accno1 << " " << two_sequences->accno2 << endl;
    #endif //DEBUG
    if (!two_sequences->aligned) {
	sequences.set_cost_matrix( 7, -5 );
	sequences.align();
    }
    #ifdef PTHREAD
    if (n_threads >1)
	pthread_mutex_lock (&databasemutex);
    #endif /* PTHREAD */
    if (two_sequences->output == 'A' || two_sequences->output == 'B') {
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
	    cerr << "Seq1 belong to cluster: " << cluster1 << ". Seq2 belong to cluster: " << cluster2 << ". Similarity: " << sequences.similarity() << endl;
	    #endif //DEBUG
	    if (sequences.similarity() > *two_sequences->cut_off) { // if (sequences.similarity(true) > *two_sequences->cut_off) {
		#ifdef DEBUG
		cerr << "Clustering seq together." << endl;
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
		    else if (cluster1.compare(cluster2)) { // No need to proceed if they are already the same
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
		    else if (cluster2.compare(cluster1)) { // No need to proceed if they are already the same
			two_sequences->db->clust_update( two_sequences->accno1, cluster2, *two_sequences->table, 1 ); // set worse sequence to point to the same sequence as better sequence
			if (!cluster1.compare( "lead" ) || !cluster1.compare( "empty" )) two_sequences->db->clust_update( two_sequences->accno1, cluster2, *two_sequences->table, 0 ); // set worse sequence to point to same sequence as better sequence
			else {
			    two_sequences->db->clust_update( cluster1, cluster2, *two_sequences->table, 1 ); // change lead of worse cluster to point to the same sequence as new best sequence
			    two_sequences->db->clust_update( cluster1, cluster2, *two_sequences->table, 0 ); // change rest of cluster to point to the same sequence as new best sequence
			}
		    }
		}
		#ifdef DEBUG
		cerr << "Clustered seq together." << endl;
		#endif //DEBUG
	    }
	    else {
		#ifdef DEBUG
		cerr << "Sequences are not clustered together." << endl;
		#endif //DEBUG
		if ( !cluster1.compare( "empty" ) ) two_sequences->db->clust_update( two_sequences->accno1, "lead", *two_sequences->table, 1 );
		if ( !cluster2.compare( "empty" ) ) two_sequences->db->clust_update( two_sequences->accno2, "lead", *two_sequences->table, 1 );
		if (!taxon_string1.empty() && taxon_string1.compare("empty") && !taxon_string2.empty() && taxon_string2.compare("empty") && (!cluster1.compare( "empty" ) || !cluster1.compare( "lead" )) && (!cluster2.compare( "empty" ) || !cluster2.compare( "lead" ))) {
		    #ifdef DEBUG
		    cerr << "Cluster 1 previous: " << cluster1 << ".  Now: " << two_sequences->db->get_cluster(two_sequences->accno1, *two_sequences->table) << endl;
		    cerr << "Cluster 2 previous: " << cluster2 << ".  Now: " << two_sequences->db->get_cluster(two_sequences->accno2, *two_sequences->table) << endl;
		    cerr << "Since sequences may get included in alignment they are counted towards the MAD: " << endl << "\t" << taxon_string1 << endl <<
			"\t" << taxon_string2 << endl << "\t" << sequences.jc_distance()-(1-sequences.similarity()) << endl;
		    #endif
		    two_sequences->deviations->insert_value( taxon_string1, taxon_string2, sequences.jc_distance()-(1-sequences.similarity()) );
		}
	    }
	}
	else {
	    #ifdef DEBUG
	    cerr << "Since sequences may get included in alignment they are counted towards the MAD: " << endl << "\t" << taxon_string1 << endl <<
		"\t" << taxon_string2 << endl << "\t" << sequences.jc_distance()-(1-sequences.similarity()) << endl;
	    #endif
	    two_sequences->deviations->insert_value( taxon_string1, taxon_string2, sequences.jc_distance()-(1-sequences.similarity()) );
	}
    }
    else if (two_sequences->outputstream != 0) {
	if ( two_sequences->output == 'a' ) {
	    if (two_sequences->output_names)
		*two_sequences->outputstream << '>' << two_sequences->accno1 << endl;
	    *two_sequences->outputstream << sequences.get_x() << endl;
	    if (two_sequences->output_names) 
		*two_sequences->outputstream << '>' << two_sequences->accno2 << endl;
	    *two_sequences->outputstream << sequences.get_y() << endl;
	}
	else if ( two_sequences->output == 'd' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) *two_sequences->outputstream << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		*two_sequences->outputstream << "Proportion sites that are different (and similarity): " << 1-sequences.similarity() << " (" << sequences.similarity() << "), Jukes-Cantor distance: " << sequences.jc_distance() << ", difference between the two: " << sequences.jc_distance()-(1.0-sequences.similarity()) << '.' << endl;
	    }
	    else {
		*two_sequences->outputstream << 1-sequences.similarity() << "/" << sequences.similarity() << "/" << sequences.jc_distance() << "/" << sequences.jc_distance()-(1.0-sequences.similarity()) << '.';
	    }
	}
	else if ( two_sequences->output == 'c' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) *two_sequences->outputstream << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		    *two_sequences->outputstream << sequences.jc_distance()-(1.0-sequences.similarity()) << endl;
	    }
	    else *two_sequences->outputstream << sequences.jc_distance()-(1.0-sequences.similarity()) << ' ';
	}
	else if ( two_sequences->output == 'j' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) *two_sequences->outputstream << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		*two_sequences->outputstream << sequences.jc_distance() << endl;
	    }
	    else *two_sequences->outputstream << sequences.jc_distance() << ' ';
	}
	else if ( two_sequences->output == 'p' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) *two_sequences->outputstream << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		*two_sequences->outputstream << 1-sequences.similarity() << endl;
	    }
	    else *two_sequences->outputstream << 1-sequences.similarity() << ' ';
	}
	else if ( two_sequences->output == 's' ) {
	    if (!two_sequences->matrix) {
		if (two_sequences->output_names) *two_sequences->outputstream << two_sequences->accno1 << " - " << two_sequences->accno2 << " | ";
		*two_sequences->outputstream << sequences.similarity() << endl;
	    }
	    else *two_sequences->outputstream << sequences.similarity() << ' ';
	}
	#ifdef DEBUG
	if (n_threads > 1) cerr << "Thread have processed: " << two_sequences->accno1 << ' ' << two_sequences->accno2 << endl;
	#endif // DEBUG
    }
    else { cerr << "No output stream!!!" << endl; }
    #ifdef PTHREAD
    if (n_threads > 1)
	pthread_mutex_unlock(&databasemutex);
    #endif /* PTHREAD */
}

void print_queue::delete_queue ( node* position ) {
    if (position != 0) {
	if (position->next != 0) delete_queue(position->next);
	delete position;
	position = 0;
    }
};

void print_queue::pop() {
    if (start != 0) {
	node* temp = start;
	start = start->next;
	delete temp;
    }
}

stringstream* print_queue::new_position () {
    if (start == 0) {
	start = new node;
	return &(start->output);
    }
    else {
	node* position = go_to_end(start);
	position->next = new node;
	return &(position->output);
    }
}

string print_queue::print_next( bool print_empty ) {
    if (start != 0) {
	#ifdef DEBUG
	cerr << "String in start: " << start->output.str() << endl;
	#endif // DEBUG
     	if (!print_empty && !something_in_first_position()) return string();
	string temp = start->output.str();
	pop();
	#ifdef DEBUG
	if (start == 0) cerr << "Everything printed!" << endl;
	cerr << "Returning string: " << temp << endl;
	#endif // DEBUG
	return temp;
    }
    else return string();
}

