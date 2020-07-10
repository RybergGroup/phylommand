/********************************************************************
Copyright (C) 2019 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

contact: martin.ryberg@ebc.uu.se
*********************************************************************/
#ifndef PARTITION_MODELS
#define PARTITION_MODELS

#include "sub_model.h"
#include "rate_model.h"
#include <vector>
#include <bitset>
#include <climits>

using namespace std;

class  partition_models {
    public:
    partition_models () { }
    partition_models ( sub_model* model ) {
	add_partition(model,NULL);
    }
    unsigned int add_partition ( sub_model* model,  rate_model* rate ) {
	return partitions.add_partition(model,rate);
    }
    unsigned int add_partition ( sub_model* model,  rate_model* rate, unsigned int first, unsigned int last ) {
	unsigned int part = add_partition ( model, rate );
	for (unsigned int i=first; i <= last; ++i) {
	    assign_pos_to_partition(i, part);
	}
    }
    unsigned int add_model_to_partition ( unsigned int partition, sub_model* model ) {
	return partitions.add_model_to_partition( partition, model);
    }
    rate_model* get_rate_model ( unsigned int pos ) {
	if (pos >= chars_to_part.size() || chars_to_part[pos] == NULL) return NULL;
	return chars_to_part[pos]->rate_mod;
    }
    sub_model* get_sub_model (unsigned int pos, unsigned int clade) {
	if (pos >= chars_to_part.size() || chars_to_part[pos] == NULL) return NULL;
	else if (clade >= chars_to_part[pos]->sub_models.size()) return NULL;
	else return chars_to_part[pos]->sub_models[0];
    }
    sub_model* get_sub_model ( unsigned int pos ) {
	get_sub_model(pos,0);
    }
    void assign_pos_to_partition ( unsigned int pos, unsigned int partition) {
	if (pos < chars_to_part.size()) {
	    chars_to_part[0] = partitions.get_item(partition);
	}
	else {
	    while (pos > chars_to_part.size()) chars_to_part.push_back(NULL);
	    chars_to_part.push_back(partitions.get_item(partition));
	}
    }
    bool check_partitions() {
	for (vector<partition_list::item*>::iterator i=chars_to_part.begin(); i != chars_to_part.end(); ++i) {
	    if (*i == 0) return false;
	    else if ((*i)->sub_models.empty()) return false;
	}
	return true;
    }
    unsigned int get_n_states(unsigned int pos) {
	return get_sub_model(pos)->get_n_states();
    }
    private:
    class partition_list {
	public:
	partition_list (): start(0) {
    
	}
	~partition_list () {
	    delete_list(start);
	}
	struct item {
	    vector<sub_model*> sub_models;
	    rate_model* rate_mod;
	    item* next;
	};
	unsigned int add_partition ( sub_model* model, rate_model* rate ) {
	    if (start == 0) {
		start = new item;
		start->sub_models.push_back(model);
		start->rate_mod = rate;
		start->next = 0;
		return 0;
	    }
	    else {
		unsigned int pos = 1;
		item* present = start;
		while (present->next != 0) {
		    present = present->next;
		    ++pos;
		}
		present->next = new item;
		present->next->sub_models.push_back(model);
		present->next->rate_mod = rate;
		present->next->next = 0;
		return pos;
	    }
	}
	unsigned int add_model_to_partition ( unsigned int partition, sub_model* model ) {
	    item* present = start;
            while (partition > 0 && present->next != 0) {
                present = present->next;
                --partition;
            }
	    if (partition > 0) cerr << "Partition search overflow" << endl;
	    present->sub_models.push_back(model);
	    return present->sub_models.size()-1;
	}
	sub_model* get_sub_model ( unsigned int partition, unsigned int clade ) {
	    item* present = start;
	    while (partition > 0 && present->next != 0) {
		present = present->next;
		--partition;
	    }
	    if (partition > 0 || clade > present->sub_models.size()) return 0;
	    else return present->sub_models[clade];
	}
	rate_model* get_rate_model ( unsigned int n ) {
	    item* present = start;
	    while (n > 0 && present->next != 0) {
		present = present->next;
		--n;
	    }
	    if (n > 0) return 0;
	    else return present->rate_mod;
	}
	item* get_item ( unsigned int partition ) {
	    item* present = start;
	    while (partition > 0 && present->next != 0) {
                present = present->next;
                --partition;
            }
	    if (partition > 0) return 0;
	    else return present;
	}
	private:
	void delete_list( item* list_item) {
	    if (list_item != 0) {
		if (list_item->next != 0) delete_list(list_item->next);
		//delete model;
		delete list_item;
	    }
	}
	item* start;
    };
    partition_list partitions;
    vector<partition_list::item*> chars_to_part;
};

#endif //PARTITION_MODELS
