#include "seqdatabase.h"

float seqdatabase::get_comp_value_pair ( const string& accno ) {
    string* sequence;
    if (accno.compare(accno1)) sequence = &sequence1;
    else if (accno.compare(accno2)) sequence = &sequence2;
    else return 0.0;
    unsigned int length = sequence->length();
    float value(0.0);
    for (unsigned int i = 0; i < length; ++i)
	if (sequence->at(i) != 'n' && sequence->at(i) != 'N') value += 1.0;
    return value;
}

void seqdatabase::move_to_next_pair_pairwisefst() {
    if (!accno1.empty()) previous_taxon = accno1;
    accno1.clear(); accno2.clear(); sequence1.clear(); sequence2.clear();
    char read_mode = '0';
    while (*input) {
	char input_char = input->get();
	if (input_char == '>') {
	    if (read_mode == '0') read_mode = 'A';
	    else if (read_mode == 'S') read_mode = 'a';
	}
	else if (input_char == '|') {
	    if (read_mode == 'A') read_mode = 'T';
	    else if (read_mode == 'a') read_mode = 't';
	}
	else if (input_char == '\n' || input_char == '\r') {
	    if (read_mode == 'A' || read_mode == 'T') read_mode = 'S';
	    else if (read_mode == 'a' || read_mode == 't') read_mode = 's';
	}
	else if (read_mode == 'A') accno1 += input_char;
	else if (read_mode == 'a') accno2 += input_char;
	else if (read_mode == 'T') taxon_strings[accno1] += input_char;
	else if (read_mode == 't') taxon_strings[accno2] += input_char;
	else if (read_mode == 'S') sequence1 += input_char;
	else if (read_mode == 's') sequence2 += input_char;
	if (read_mode == 's' && (input->peek() == '>' || input->peek() == EOF)) break;
    }
    #ifdef DEBUG
    cerr << previous_taxon << " : " << accno1 << endl;
    #endif //DEBUG
    if (!previous_taxon.compare("empty") || previous_taxon.compare(accno1)) mode = '2';
    else mode = '3';
    if (input->bad() || input->peek() == EOF) mode = '9';
}

void seqdatabase::move_to_next_pair_fst (bool only_lead) {
    #ifdef DEBUG
    cerr << "Moving to next pair" << endl;
    #endif //DEBUG
    if (mode != '9') { // unless we should quit
	while (1) {
	    if (mode == '0') { // if at start of new round
		if (previous_taxon.compare("empty")) { // if no at first run
		    if (!fst.set_seq1(previous_taxon)) mode = '9'; // if we fail to set first taxon to previous taxon set to quit
		    if (!fst.next_seq1()) mode = '9'; // if not able to go to next sequence set to quit
		}
		else { // if at first round
		    if (!fst.initiate_sequence_retrieval()) mode = '9'; // set to quit is we fail to initiate
		}
		if (mode != 9) {
		    fst.set_seq2_to_seq1(); //set seq 2 to to be same as seq 1 i.e. previous taxon
		    accno2.clear();
		    sequence2.clear();
		    accno1 = fst.get_accno1(); // get accno
		    fst.get_sequence1(sequence1); // get sequence
		    previous_taxon = accno1; // set previous accno to be present accno
		    mode = '1'; // we have now read seq1
		    #ifdef DEBUG
		    cerr << "New first seq: " << accno1 << endl;
		    #endif //DEBUG
		}
	    }
	    else if (mode != '9') { // if not at end
		if (fst.next_seq2()) { // if able to read next sequence
		    accno2 = fst.get_accno2(); // get accno
		    fst.get_sequence2(sequence2); // get sequence
		    if (mode == '1') {
			if (fst.seq2_is_last_seq()) mode = '9'; // if at last seq quit
			else mode = '2'; // we have now read a second sequence
		    }
		    else if (mode == '2') mode = '3';
		    #ifdef DEBUG
		    cerr << "New second seq: " << accno2 << endl;
		    #endif //DEBUG
		    break;
		}
		else if (mode == '1') { mode = '9'; break; }// if not able to read a second sequence we are at the end
		else mode = '0'; // if second sequence read but there are no more, start a new run
	    }
	}
    }
}

#ifdef DATABASE
bool seqdatabase::alignment_groups_present_sql() {
    if (OPEN) {
	sqlite3_stmt* statement;
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

bool seqdatabase::create_alignment_groups_sql(){
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

vector<string> seqdatabase::tables_in_database_sql() {
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

bool seqdatabase::initiate_sequence_retrieval_sql ( const string table, string min_length) {
    mode = '0';
    query = "SELECT accno,sequence,cluster FROM ";
    query += table;
    query += " WHERE LENGTH(sequence)>";
    query += min_length;
    query += ';';
    if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) == SQLITE_OK) return true;
    else return false;
}

void seqdatabase::move_to_next_pair_sql (bool only_lead) {
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
	    mode = '1';
	}
	else mode = '9';
    }
    else if (mode == '0') mode = '1';
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
		if (mode == '2' || mode == '1' || mode == '0') mode = '9'; // if on first or second sequence, all pairs have been done
		else if(sqlite3_prepare_v2(db, query.c_str(), -1, &statement, 0) != SQLITE_OK) mode = '9'; // if failing sql statment
		else mode = 0; // otherwise we start from the begining again
		break;
	    }
	}
    }
}


bool seqdatabase::insert_alignment_group_sql (const string& table, const string& group) {
    if (OPEN) {
	if (group.compare("empty")) { // if something in alignment groups
	    //sqlite3_stmt* statement;
    	    string taxon;
	    char alignable;
	    for (unsigned int i=0; i < group.length(); ++i) {
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
	return true;
    }
    else return false;
}

void seqdatabase::get_taxon_string_sql( string accno, string& taxon_string ) {
    sqlite3_stmt *statement;
    //string taxon_string;
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
    //return taxon_string;
}

void seqdatabase::clust_update_sql( const string accno, const string cluster, const string table, const bool where_accno) {
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

string seqdatabase::get_cluster_sql( const string accno, const string table ) {
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

float seqdatabase::get_comp_value_sql ( const string& accno, const string& table ) {
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
