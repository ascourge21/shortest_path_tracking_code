/*
 * matrix3d.h
 *
 *  Created on: Sep 9, 2015
 *      Author: nripesh
 */

#ifndef MATRIX3D_H_
#define MATRIX3D_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>

using namespace boost::numeric::ublas;

class matrix3d {
public:
	matrix3d();
	matrix3d(int N, int M, int L); // matrix of size: rows N, cols M, depth L, with 0 default
	void assign(int i, int j, int k, double val);
	double at(int i, int j, int k);
	std::vector<double> at(int i, int j); // get vector at row (i, j)
	matrix<double> atRow(int row); // get mat at row i 
	int rowSize();
	int colSize();
	int depth();
	void print();
	virtual ~matrix3d();
private:
	int N;
	int M;
	int L;
	matrix<std::vector<double> > mat;
	std::vector<matrix<double> > matByRow;
};

#endif /* MATRIX3D_H_ */
