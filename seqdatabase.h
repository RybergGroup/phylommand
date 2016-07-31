#include <iostream>
#include <stdlib.h>
#ifdef DATABASE
#include "sqlite3.h"
#endif //DATABASE
#include <string.h>
#include <vector>
#include <indexedfasta.h>

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
    #ifdef DATABASE
    void open_sqldatabase( const char* name ) {
	if (sqlite3_open(name, &db) == 0 ) OPEN = true;
	else OPEN = false;
    };
    bool alignment_groups_present ();
    bool create_alignment_groups();
    vector<string> tables_in_database();
    bool insert_alignment_group (const string& table, const string& group);
    void clust_update( const string accno, const string cluster, const string table, const bool where_accno);
    string get_cluster( const string accno, const string table );
    float get_comp_value ( const string accno, const string table);
    string get_taxon_string( string accno );
    bool initiate_sequence_retrieval( const string table, string min_length);
    void move_to_next_pair() { move_to_next_pair(false); };
    void move_to_next_pair( bool only_lead );
    #endif //DATABASE
    //// Corresponding functions for fasta
    void open_fastafile( const char* name ) {
	fst.open(name);
	if (fst.is_open()) OPEN = true;
	else OPEN = false;
    };
    bool alignment_groups_present_fst() { return alignment_groups.open(); };
    bool create_alignment_groups() {
	string name = databasename + "alignment_groups";
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
    string get_taxon_string_fst ( const string accno ) { fst.get_taxon_string( accno ); };
    bool initiate_sequence_retrival_fst ( return fst.index(); };
    void move_to_next_pair_fst ( bool only_lead ) {
	// needs to be written
	// will depend probably depend on several indexedfasta functions
	// should assign values to accno1 and 2, and sequence1 and 2
    };
    //////////////////////////////////////////
    bool all_pairs() { if(mode == '9') return true; else return false; };
    string get_accno1() { return accno1; };
    string get_accno2() { return accno2; };
    string get_sequence1() { return sequence1; };
    string get_sequence2() { return sequence2; };
////////////////////////////////////////////////////////
    private:
    class cluster_struct {
	public:
	void update ( const string accno, const string cluster, const bool where_accno); // needs to be written
	string get_cluster( const string accno ); // needs to be written
	private:
    };
    string databasename;
    char databasetype;
    bool OPEN;
    char mode;
    //////// For fasta files
    indexedfasta fst;
    ofstream alignment_groups;
    cluster_struct clusters;
    /////////////////////
    string previous_taxon;
    string accno1;
    string sequence1;
    string accno2;
    string sequence2;
    #ifdef DATABASE
    sqlite3 *db;
    sqlite3_stmt *statement;
    string query;
    #endif //DATABASE
};
