#include <iostream>
#include <vector>
#include "set.h"

set::set(){
	list = std::vector<int>();
}

void set::insert(int a) {
	// if (!(std::find(list.begin(), list.end(), a) != list.end())) { // if exists
	// 	list.push_back(a);
	// }
	bool found = false;
	for (int i = 0; i < list.size(); i++) {
		if (list[i] == a) found = true;
	}
	if (!found) list.push_back(a);
}

void set::insert(std::vector<int> v) {
	for (int i = 0; i < v.size(); i++) {
		insert(v[i]);
	}
}

void set::insert(std::vector<double> v) {
	for (int i = 0; i < v.size(); i++) {
		insert((int)v[i]);
	}
}

int set::get_count() {
	return list.size();
}

std::vector<int> set::get_elements() {
	return list;
}
