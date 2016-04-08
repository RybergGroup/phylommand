//#include <stdlib.h>
#include <iostream>
#include "sqlite3.h"

using namespace std;

main (int argc, char *argv[]) {
    sqlite3 *db;
    cout << sqlite3_open("../inocybe/inocyb.db", &db) << endl;
    sqlite3_stmt *statement;
    cout << sqlite3_prepare_v2(db, "SELECT species FROM gb_data WHERE key='JN642244'", -1, &statement, 0) << endl;
    
    cout << sqlite3_close(db) << endl;
}
