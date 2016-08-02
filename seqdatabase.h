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
    seqdatabase(): databasetype('f'), OPEN(false), previous_taxon("empty"), mode(0) {};
    seqdatabase(const char* file): databasetype('f'), previous_taxon("empty"), mode(0), databasename(file) {
    	databasetype = 'f';
	open_fastafile( file );
    }
    seqdatabase(const char* file, const char* type): previous_taxon("empty"), mode(0), databasename(file) {
	if (!strcmp(type,"fasta")){
	    databasetype = 'f';
	    open_fastafile( file );
	}
	#ifdef DATABASE
	else if (!strcmp(type,"sqlite")) {
	    databasetype = 's';
	    open_sqldatabase( file );
	}
	#endif //DATABASE
	else {
	    databasetype = 'e';
	    OPEN = false;
	}
    }
    bool alignment_groups_present () {
	if (databasetype == 'f') return alignment_groups_present_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return alignment_groups_present_sql ();
	#endif //DATABASE
	else return false;
    };
    bool create_alignment_groups() {
	if (databasetype == 'f') return create_alignment_groups_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return create_alignment_groups_sql ();
	#endif //DATABASE
	else return false;
    };
    vector<string> tables_in_database() {
	if (databasetype == 'f') return tables_in_database_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return tables_in_database_sql ();
	#endif //DATABASE
	else {
	    vector<string> temp; return temp;
	}
    };
    bool insert_alignment_group (const string& table, const string& group) {
	if (databasetype == 'f') return insert_alignment_group_fst (table, group);
	#ifdef DATABASE
	else if (databasetype == 's') return insert_alignment_group_sql (table, group);
	#endif //DATABASE
	else return false;
    };
    void clust_update( const string accno, const string cluster, const string table, const bool where_accno) {
	if (databasetype == 'f') return clust_update_fst (accno, cluster, table, where_accno);
	#ifdef DATABASE
	else if (databasetype == 's') return clust_update_sql (accno, cluster, table, where_accno);
	#endif //DATABASE
    };
    string get_cluster ( const string accno, const string table ) {
	if (databasetype == 'f') return get_cluster_fst (accno, table);
	#ifdef DATABASE
	else if (databasetype == 's') return get_cluster_sql (accno, table);
	#endif //DATABASE
	else { string temp; return temp; }
    };
    float get_comp_value ( const string accno, const string table) {
	if (databasetype == 'f') return get_comp_value_fst (accno, table);
	#ifdef DATABASE
	else if (databasetype == 's') return get_comp_value_sql ( accno, table);
	#endif //DATABASE
	else return 0.0;
    };
    string get_taxon_string( string accno ) {
	string return_string;
	if (databasetype == 'f') get_taxon_string_fst ( accno, return_string );
	#ifdef DATABASE
	else if (databasetype == 's') get_taxon_string_sql ( accno, return_string );
	#endif //DATABASE
	return return_string;
    };
    bool initiate_sequence_retrieval( const string table, string min_length) {
	if (databasetype == 'f') return initiate_sequence_retrieval_fst ();
	#ifdef DATABASE
	else if (databasetype == 's') return initiate_sequence_retrieval_sql ( table, min_length);
	#endif //DATABASE
	else return false;
    };
    void move_to_next_pair( bool only_lead ) {
	if (databasetype == 'f') move_to_next_pair_fst ( only_lead );
	#ifdef DATABASE
	else if (databasetype == 's') move_to_next_pair_sql ( only_lead );
	#endif //DATABASE
    };
    bool all_pairs() { if(mode == '9') return true; else return false; };
    string get_accno1() { return accno1; };
    string get_accno2() { return accno2; };
    string get_sequence1() { return sequence1; };
    string get_sequence2() { return sequence2; };
    void print_clusters() {
	string name = databasename + ".clusters";
	ofstream output (name.c_str(), ios::out);
	clusters.print_clusters( output );
	output.close();
    }
////////////////////////////////////////////////////////
    private:
    class single_link_clusters {
	public:
	void update ( const string accno, const string cluster, const bool where_accno) {
	    if (where_accno) {
		if (!cluster.compare("lead")) clusters[accno] = set<string>();
		else if (clusters.find(cluster) != clusters.end()) {
		    clusters[cluster].insert(accno);
		}
	    }
	    else {
		if (clusters.find(accno) != clusters.end()) {
		    if (clusters.find(cluster) != clusters.end()) {
			clusters[cluster].insert(clusters[accno].begin(),clusters[accno].end());
			clusters.erase(accno);
		    }
		    else {
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
    /////////////////////
    string previous_taxon;
    string accno1;
    string sequence1;
    string accno2;
    string sequence2;
    //// functions for fasta
    void open_fastafile( const char* name ) {
	fst.open(name);
	if (fst.is_open()) OPEN = true;
	else OPEN = false;
    };
    bool alignment_groups_present_fst() { return alignment_groups.is_open(); };
    bool create_alignment_groups_fst() {
	string name = databasename + ".alignment_groups";
	alignment_groups.open(name.c_str());
	return alignment_groups.is_open();
    };
    vector<string> tables_in_database_fst() {
	vector<string> temp;
	temp.push_back(databasename);
	return temp;
    };
    bool insert_alignment_group_fst (const string& table, const string& group) { alignment_groups << table << "\t" << group << endl; };
    void clust_update_fst ( const string accno, const string cluster, const string table, const bool where_accno) {
	clusters.update( accno, cluster, where_accno );
    };
    string get_cluster_fst( const string accno, const string table ) {
	return clusters.get_cluster( accno );
    };
    float get_comp_value_fst ( const string accno, const string table) { fst.get_comp_value( accno ); }; // needs to be written
    void get_taxon_string_fst ( const string accno, string& return_string ) { fst.get_taxon_string( accno, return_string ); };
    bool initiate_sequence_retrieval_fst () { mode = '0'; return fst.initiate_sequence_retrieval(); };
    void move_to_next_pair_fst ( bool only_lead );
    //////////////////////////////////////////
    #ifdef DATABASE
    sqlite3 *db;
    sqlite3_stmt *statement;
    string query;
    void open_sqldatabase( const char* name ) {
	if (sqlite3_open(name, &db) == 0 ) OPEN = true;
	else OPEN = false;
    };
    bool alignment_groups_present_sql ();
    bool create_alignment_groups_sql();
    vector<string> tables_in_database_sql();
    bool insert_alignment_group_sql (const string& table, const string& group);
    void clust_update_sql( const string accno, const string cluster, const string table, const bool where_accno);
    string get_cluster_sql( const string accno, const string table );
    float get_comp_value_sql ( const string accno, const string table);
    void get_taxon_string_sql ( string accno, string& taxon_string );
    bool initiate_sequence_retrieval_sql( const string table, string min_length);
    void move_to_next_pair_sql () { move_to_next_pair(false); };
    void move_to_next_pair_sql ( bool only_lead );
    #endif //DATABASE
};
