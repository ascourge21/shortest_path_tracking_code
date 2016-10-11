/*
 * shortestpaths_3d.h
 *
 *  Created on: Jan 2016
 *      Author: nripesh
 */

#ifndef SHORTEST_PATHS_3D_FAST_H_
#define SHORTEST_PATHS_3D_FAST_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "matrix3d.h"
#include "math.h"
#include "mex.h"
#include "set.h"

using namespace boost::numeric::ublas;

class shortest_paths_3d_fast {
public:
	shortest_paths_3d_fast(int fr_total = 40, int fr_es = 17, std::vector<double> = NULL);
    
    // the one that only uses regular norm -> 2nd order included
	void findShortestPaths(matrix<double> &x, matrix<double> y, matrix<double> &z, 
			matrix<double> &init_xyz, matrix3d &near_nbor_inds,  
			matrix3d &nbors_dist_Mat, matrix3d &nbors_shp_Mat,
			std::vector<int> &inds_to_track, std::vector<int> &init_ref_inds, 
            matrix<int> &path_inds, std::vector<double> &path_dist);	// the function that finds shortest path findShortestPaths();

	// destructor
	virtual ~shortest_paths_3d_fast();
private:

	int fr_total; 
	int fr_es; 
	std::vector<double> fr_wts; 
    
    double DIST_IN_WT, DIST_PREV_WT, PT_DERIV_WT, SH_PREV_WT ;

	// vector operations
	double norm(std::vector<double> &a, std::vector<double> &b); // distance between a and b -> norm
    std::vector<double> subtract(std::vector<double> &a, std::vector<double> &b);
    std::vector<double> add(std::vector<double> &a, std::vector<double> &b);
    std::vector<double> multiply(std::vector<double> &a, double b);
   
    double dist_new(std::vector<double> &current_pt, std::vector<double> &root_pt,
			std::vector<double> &initial_pt, std::vector<double> &root_root_pt, 
			double root_to_pt_euc_dist, double root_to_pt_shp_dist,
            int fr_i);
   
	int minPosInRow(matrix<double> &m, int col);
};

#endif /* SHORTEST_PATHS_3D_FAST_H_ */
