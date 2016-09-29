#include <iostream>
#include <stdlib.h>
#include <map>
#include <set>
#ifdef DATABASE
#include "sqlite3.h"
#endif //DATABASE
#include <string.h>
#include <vector>
#include "indexedfasta.h"

using namespace std;

class seqdatabase {
    public:
    seqdatabase(): databasetype('f'), OPEN(false), mode(0), previous_taxon("empty") {};
    seqdatabase(const string file): databasename(file), databasetype('f'), mode(0), previous_taxon("empty") {
    	databasetype = 'f';
	open_fastafile( );
    }
    seqdatabase(const string file, const char* type): databasename(file), mode(0), previous_taxon("empty") {
	if (!strcmp(type,"fasta")){
	    databasetype = 'f';
	    open_fastafile( );
	}
	else if (!strcmp(type,"pairfa")) {
	    databasetype = 'P';
	    open_file( );
	}
	#ifdef DATABASE
	else if (!strcmp(type,"sqlite")) {
	    databasetype = 's';
	    open_sqldatabase( );
	}
	#endif //DATABASE
	else {
	    databasetype = 'e';
	    OPEN = false;
	}
    }
    bool alignment_groups_present () {
	if (databasetype == 'f' || databasetype == 'P') return alignment_groups_present_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return alignment_groups_present_sql ();
	#endif //DATABASE
	else return false;
    };
    bool create_alignment_groups() {
	if (databasetype == 'f' || databasetype == 'P') return create_alignment_groups_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return create_alignment_groups_sql ();
	#endif //DATABASE
	else return false;
    };
    vector<string> tables_in_database() {
	if (databasetype == 'f' || databasetype == 'P') return tables_in_database_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return tables_in_database_sql ();
	#endif //DATABASE
	else {
	    vector<string> temp; return temp;
	}
    };
    bool insert_alignment_group (const string& table, const string& group) {
	if (databasetype == 'f' || databasetype == 'P') return insert_alignment_group_fst (table, group);
	#ifdef DATABASE
	else if (databasetype == 's') return insert_alignment_group_sql (table, group);
	#endif //DATABASE
	else return false;
    };
    void clust_update( const string accno, const string cluster, const string table, const bool where_accno) {
	if (databasetype == 'f' || databasetype == 'P') return clust_update_fst (accno, cluster, table, where_accno);
	#ifdef DATABASE
	else if (databasetype == 's') return clust_update_sql (accno, cluster, table, where_accno);
	#endif //DATABASE
    };
    string get_cluster ( const string accno, const string table ) {
	if (databasetype == 'f' || databasetype == 'P') return get_cluster_fst (accno, table);
	#ifdef DATABASE
	else if (databasetype == 's') return get_cluster_sql (accno, table);
	#endif //DATABASE
	else { string temp; return temp; }
    };
    float get_comp_value ( const string& accno, const string& table) {
	if (databasetype == 'f') return get_comp_value_fst (accno, table);
	else if (databasetype == 'P') return get_comp_value_pair (accno);
	#ifdef DATABASE
	else if (databasetype == 's') return get_comp_value_sql ( accno, table);
	#endif //DATABASE
	else return 0.0;
    };
    string get_taxon_string( string accno ) {
	string return_string;
	if (databasetype == 'f') get_taxon_string_fst ( accno, return_string );
	if (databasetype == 'P') get_taxon_string_from_map( accno, return_string);
	#ifdef DATABASE
	else if (databasetype == 's') get_taxon_string_sql ( accno, return_string );
	#endif //DATABASE
	return return_string;
    };
    bool initiate_sequence_retrieval( const string table, string min_length) {
	if (databasetype == 'f') return initiate_sequence_retrieval_fst ();
	if (databasetype == 'P') return initiate_sequence_retrival_file();
	#ifdef DATABASE
	else if (databasetype == 's') return initiate_sequence_retrieval_sql ( table, min_length);
	#endif //DATABASE
	else return false;
    };
    void move_to_next_pair( bool only_lead ) {
	if (databasetype == 'f') move_to_next_pair_fst ( only_lead );
	if (databasetype == 'P') move_to_next_pair_pairwisefst ( );
	#ifdef DATABASE
	else if (databasetype == 's') move_to_next_pair_sql ( only_lead );
	#endif //DATABASE
    };
    bool at_new_first() { if (mode == '0' || mode == '1' || mode == '2' || mode == '9') return true; else return false; }
    //bool end_of_round() { if (mode == '0' || mode == '9') return true; else return false; }
    bool all_pairs() { if(mode == '9') return true; else return false; };
    string get_accno1() { return accno1; };
    string get_accno2() { return accno2; };
    string get_sequence1() { return sequence1; };
    string get_sequence2() { return sequence2; };
    void print_clusters(ostream& output) {
	//string name; 
	//if (!databasename.empty()) name = databasename + ".clusters";
	//else name = "sequence.clusters";
	//ofstream output (name.c_str(), ios::out);
	clusters.print_clusters( output );
	//output.close();
    }
    void add_taxonomy( map<string,string>& taxonomy ) {
	taxon_strings = taxonomy;
    };
////////////////////////////////////////////////////////
    private:
    class single_link_clusters {
	public:
	void update ( const string accno, const string cluster, const bool where_accno) {
	    if (!accno.compare(cluster)) {
		cerr << "Warning!!! Accession number cluster name collision when updating clusters." << endl;
		return;
	    }
	    #ifdef DEBUG
	    cerr << "Adding " << accno << " to " << cluster << endl;
	    #endif //DEBUG
	    if (where_accno) {
		#ifdef DEBUG
		if (!cluster.compare("lead")) cerr << "Making cluster for " << accno << endl;
		else if (clusters.find(cluster) != clusters.end()) cerr << "Adding " << accno << " to cluster for " << cluster << endl; 
		#endif //DEBUG
		if (!cluster.compare("lead")) clusters[accno] = set<string>();
		else if (clusters.find(cluster) != clusters.end()) {
		    clusters[cluster].insert(accno);
		}
	    }
	    else { // Where cluster
		if (clusters.find(accno) != clusters.end()) {
		    #ifdef DEBUG
		    cerr << "Cluster " << accno << " exist" << endl;
		    #endif //DEBUG
		    if (clusters.find(cluster) != clusters.end()) {
			#ifdef DEBUG
			cerr << "\tAdding cluster " << accno << " to cluster " << cluster<< endl;
			#endif //DEBUG
			clusters[cluster].insert(clusters[accno].begin(),clusters[accno].end());
			clusters.erase(accno);
		    }
		    else {
			#ifdef DEBUG
			cerr << "\tNo cluster for " << cluster << ", creating it and adding cluster " << accno << endl;
			#endif //DEBUG
			clusters[cluster] = clusters[accno];
			clusters.erase(accno);
		    }
		}
    	    }
	}; //
	string get_cluster( const string accno ) {
	    string return_string = "lead";
	    if (clusters.find(accno) != clusters.end()) return return_string;
	    for (map<string, set<string> >::iterator i=clusters.begin(); i != clusters.end(); ++i) {
		if (i->second.find(accno) != i->second.end()) return i->first;
	    }
	    return_string = "empty";
	    return return_string;
	}; //
	void print_clusters (ostream& output) {
	    for (map<string, set<string> >::iterator i=clusters.begin(); i != clusters.end(); ++i) {
		output << i->first;
		for (set<string>::iterator j= i->second.begin(); j!= i->second.end(); ++j) output << ' ' << *j;
		output << endl;
	    }
	}
	private:
	map<string,set<string> > clusters; // name sequences
    };
    string databasename;
    char databasetype;
    bool OPEN;
    char mode;
    //////// For fasta files
    indexedfasta fst;
    ofstream alignment_groups;
    single_link_clusters clusters;
    //////// For pairwise fasta
    ifstream file;
    istream* input;
    map<string,string> taxon_strings;
    /////////////////////
    string previous_taxon;
    string accno1;
    string sequence1;
    string accno2;
    string sequence2;
    //// functions for fasta
    void open_fastafile( ) {
	fst.open(databasename);
	if (fst.is_open()) OPEN = true;
	else OPEN = false;
    };
    bool alignment_groups_present_fst() { return alignment_groups.is_open(); };
    bool create_alignment_groups_fst() {
	string name; 
	if (!databasename.empty()) name = databasename + ".alignment_groups";
	else name = "sequence.alignment_groups";
	alignment_groups.open(name.c_str());
	return alignment_groups.is_open();
    };
    vector<string> tables_in_database_fst() {
	vector<string> temp;
	temp.push_back(databasename);
	return temp;
    };
    bool insert_alignment_group_fst (const string& table, const string& group) {
	if (alignment_groups.good()) {
	    alignment_groups << table << "\t" << group << endl;
	    return true;
	}
	else return false;
    };
    void clust_update_fst ( const string accno, const string cluster, const string table, const bool where_accno) {
	clusters.update( accno, cluster, where_accno );
    };
    string get_cluster_fst( const string accno, const string table ) {
	return clusters.get_cluster( accno );
    };
    float get_comp_value_fst ( const string& accno, const string& table) { return fst.get_comp_value( accno ); }; 
    void get_taxon_string_fst ( const string accno, string& return_string ) {
	if (!taxon_strings.empty()) get_taxon_string_from_map(accno, return_string);
	else fst.get_taxon_string( accno, return_string );
    };
    bool initiate_sequence_retrieval_fst () { mode = '0'; return fst.initiate_sequence_retrieval(); };
    void move_to_next_pair_fst ( bool only_lead );
    //////////////////////////////////////////
    //////// Funktions for pairwise fasta
    void open_file () {
	if (databasename.empty()) {
	    input = &cin;
	    OPEN = true;
	}
	else {
	    file.open(databasename.c_str());
	    if (file.good()) {
		OPEN = true;
		input = &file;
	    }
	}
    }
    float get_comp_value_pair ( const string& accno );
    void get_taxon_string_from_map ( const string& accno, string& return_string) {
	map<string,string>::iterator i = taxon_strings.find(accno);
	if (i != taxon_strings.end()) return_string = i->second;
	else {
	    i = taxon_strings.find("default");
	    if (i != taxon_strings.end()) return_string = i->second;
	}
    }
    bool initiate_sequence_retrival_file() {
	mode = 0;
	if ( input != &cin) { // alternative _sCheck(input != &cin)
	    if (file.bad()) file.clear();
	    file.seekg(0,file.beg);
	    return file.good();
	}
	else return true;
    }
    void move_to_next_pair_pairwisefst();
    #ifdef DATABASE
    sqlite3 *db;
    sqlite3_stmt *statement;
    string query;
    void open_sqldatabase( ) {
	if (sqlite3_open(databasename.c_str(), &db) == 0 ) OPEN = true;
	else OPEN = false;
    };
    bool alignment_groups_present_sql ();
    bool create_alignment_groups_sql();
    vector<string> tables_in_database_sql();
    bool insert_alignment_group_sql (const string& table, const string& group);
    void clust_update_sql( const string accno, const string cluster, const string table, const bool where_accno);
    string get_cluster_sql( const string accno, const string table );
    float get_comp_value_sql ( const string& accno, const string& table);
    void get_taxon_string_sql ( string accno, string& taxon_string );
    bool initiate_sequence_retrieval_sql( const string table, string min_length);
    void move_to_next_pair_sql () { move_to_next_pair_sql(false); };
    void move_to_next_pair_sql ( bool only_lead );
    #endif //DATABASE
};
