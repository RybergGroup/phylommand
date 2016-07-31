#include "seqdatabase.h"

#ifdef DATABASE
bool seqdatabase::alignment_groups_present() {
    if (OPEN) {
	sqlite3_stmt* statment;
	const char query[] = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;";
	bool flag = 0;
	if(sqlite3_prepare_v2(db, query, -1, &statement, 0) == SQLITE_OK) {
	    while (sqlite3_step(statement) == SQLITE_ROW) {
		string table = (char*)sqlite3_column_text(statement,0);
		if (!table.compare("alignment_groups")) {
		    flag = 1;
		    break;
		}
	    }
	}
	sqlite3_finalize(statement);
	return flag;
    }
    else return false;
}

bool seqdatabase::create_alignment_groups(){
    if (OPEN) {
        sqlite3_stmt *statement;
        const char query[] = "CREATE TABLE alignment_groups (gene TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', tree TEXT DEFAULT 'empty', tree_method TEXT DEFAULT 'empty', alignable INTEGER, PRIMARY KEY (gene, taxon));";
        if(sqlite3_prepare_v2(db, query, -1, &statement, 0) == SQLITE_OK) {
            if (sqlite3_step(statement) != SQLITE_DONE) {
                std::cerr << "Could not create table for alignment_groups (1)." << endl;
		sqlite3_finalize(statement);
                return false;
            }
            sqlite3_finalize(statement);
	    return true;
        }
	else return false;
    }
    else {
	cerr << "No connection to database. Can not create alignment_group table." << endl;
	return false;
    }
}

vector<string> seqdatabase::tables_in_database() {
    const char query[] = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;";
    sqlite3_stmt *statement;
    vector<string> return_vector;
    if (OPEN) {
	if (sqlite3_prepare_v2(db, query, -1, &statement, 0) == SQLITE_OK) {
	    while (sqlite3_step(statement) == SQLITE_ROW) {
		return_vector.push_back( (char*)sqlite3_column_text(statement,0) );
	    }
	}
	sqlite3_finalize(statement);
    }
    return return_vector;
}

bool seqdatabase::initiate_sequence_retrieval( const string table, string min_length) {
    query = "SELECT accno,sequence,cluster FROM ";
    query += table;
    query += " WHERE LENGTH(sequence)>";
    query += min_length;
    query += ';';
    if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) return true;
    else return false;
}

void seqdatabase::move_to_next_pair(bool only_lead) {
    if (mode=='0' && previous_taxon.compare("empty")) {
	if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
	    while (1) { // Move up until previou sequence, we are finished will all before that and should start from next
		string accno;
		if (sqlite3_step(statement) == SQLITE_ROW) accno = (char*)sqlite3_column_text(statement,0);
		else {
		    mode = '9';
		    break;
		}
		if (!accno.compare(previous_taxon)) break;
	    }
	    mode = 1;
	}
	else mode = '9';
    }
    if (mode != '9') {
	while (1) {
	    if(sqlite3_step(statement) == SQLITE_ROW) {
		if (only_lead) { // if only doing lead sequences
		    string cluster = (char*)sqlite3_column_text(statement,2);
		    if (cluster.compare("lead")) { // if not lead sequence move on
			if (mode == '1') previous_taxon = (char*)sqlite3_column_text(statement,0); // if looking for first seq save accno
			continue; // if not lead sequence continue
		    }
		}
		if (mode == '1') {
		    accno1 = (char*)sqlite3_column_text(statement,0);
		    previous_taxon = accno1;
		    sequence1 = (char*)sqlite3_column_text(statement,1);
		    mode = '2';
		}
		else if ( mode == '2' || mode == '3') {
		    accno2 = (char*)sqlite3_column_text(statement,0);
		    sequence2 = (char*)sqlite3_column_text(statement,0);
		    if (mode == '2') mode = '3';
		    break;
		}
	    }
	    else { // reached the end
		sqlite3_finalize(statement);
		if (mode == '2' || mode == '1') mode = '9'; // if on first or second sequence, all pairs have been done
		else if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) != SQLITE_OK) mode = '9'; // if failing sql statment
		else mode = 0; // otherwise we start from the begining again
		break;
	    }
	}
    }
}

bool seqdatabase::insert_alignment_group (const string& table, const string& group) {
    if (OPEN) {
	if (group.compare("empty")) { // if something in alignment groups
	    sqlite3_stmt* statment;
    	    string taxon;
	    char alignable;
	    for (int i=0; i < group.length(); ++i) {
		if (group[i] == ';') {
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
		else if (group[i] == '_') {
		    if (group[++i] == 'A') alignable = '1';
		    else alignable = '0';
		}
		else taxon += group[i];
	    }
	}
    }
    else return false;
}

string seqdatabase::get_taxon_string( string accno ) {
    sqlite3_stmt *statement;
    string taxon_string;
    if (OPEN) {
	string query = "SELECT taxon_string FROM ";
	query += "gb_data";
	query += " WHERE accno='";
	query += accno;
	query += "';";
	if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
	    if (sqlite3_step(statement) == SQLITE_ROW) {
		taxon_string= (char*)sqlite3_column_text(statement,0);
	    }
	}
	sqlite3_finalize(statement);
    }
    return taxon_string;
}

void seqdatabase::clust_update( const string accno, const string cluster, const string table, const bool where_accno) {
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

string seqdatabase::get_cluster( const string accno, const string table ) {
    string cluster;
    if (OPEN) {
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
    }
    return cluster;
}

float seqdatabase::get_comp_value ( const string accno, const string table ) {
    int length=0;
    float prop_N=0.0;
    if (OPEN) {
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
    }
    return length*(1-prop_N);
}

#endif //DATABASE
