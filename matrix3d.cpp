/*
 * matrix3d.cpp
 *
 *  Created on: Sep 9, 2015
 *      Author: nripesh
 */

#include "matrix3d.h"

// default constructor
matrix3d::matrix3d() {
	N = 0;
	M = 0;
	L = 0;
	mat = matrix<std::vector<double> > ();
	matByRow = std::vector<matrix<double> > ();
}

// mat should have been initialized before
matrix3d::matrix3d(int N, int M, int L) {
	// TODO Auto-generated constructor stub

	mat = matrix<std::vector<double> >(N, M);
	matByRow = std::vector<matrix<double> > (N);
	for (int i = 0; i < N; i++) {
		matByRow[i] = matrix<double>(M, L);
		for (int j = 0; j < M; j++) {
			mat(i, j) = std::vector<double>(L);
		}
	}

	this->N = N;
	this->M = M;
	this->L = L;
}

void matrix3d::assign(int i, int j, int k, double val) {
	if (i < 0 || i > N-1)
		throw std::out_of_range("assign: i out of range");
	if (j < 0 || j > M-1)
		throw std::out_of_range("assign: j out of range");
	if (k < 0 || k > L-1)
		throw std::out_of_range("assign: k out of range");
	mat(i, j)[k] = val;
	matByRow[i](j, k) = val;
}

double matrix3d::at(int i, int j, int k) {
	if (i < 0 || i > N-1)
		throw std::out_of_range("at: i out of range");
	if (j < 0 || j > M-1)
		throw std::out_of_range("at: j out of range");
	if (k < 0 || k > L-1)
		throw std::out_of_range("at: k out of range");
	return mat(i, j)[k];
}

std::vector<double> matrix3d::at(int i, int j) {
	if (i < 0 || i > N-1)
		throw std::out_of_range("at vec: i out of range");
	if (j < 0 || j > M-1)
		throw std::out_of_range("at vec: j out of range");
	return mat(i, j);
}

matrix<double> matrix3d::atRow(int row) {
	if (row < 0 || row > N-1)
		throw std::out_of_range("at row: row out of range");
	return matByRow[row];
}


void matrix3d::print() {
	for (int k = 0; k < L; k++) {
		std::cout << "level: "<< k << std::endl;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				std::cout << mat(i, j)[k] << ", ";
			}
			std::cout << std::endl;
		}
	}
}


int matrix3d::rowSize() {
	return N;
}

int matrix3d::colSize() {
	return M;
}

int matrix3d::depth() {
	return L;
}


matrix3d::~matrix3d() {
	// TODO Auto-generated destructor stub
}

