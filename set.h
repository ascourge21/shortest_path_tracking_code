#include <iostream>
#include <vector>
#ifndef SET_H_
#define SET_H_

class set {
	public:
		set();
		void insert(int a);
		void insert(std::vector<int> v);
		void insert(std::vector<double> v);	
		int get_count();
		virtual ~set(){};
		std::vector<int> get_elements();

	private:
		std::vector<int> list;
		int count;
};

#endif /* SET_H_ */