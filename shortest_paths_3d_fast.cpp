/*
 * shortest_paths_3d_fast.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: nripesh
 */



#include "shortest_paths_3d_fast.h"
#include <stdlib.h>
#include <sstream>

// to convert string to double, because c++'s normal thing didn't work with MEX.
double string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
	 return 0;
	return x;
} 


/**
 */
shortest_paths_3d_fast::shortest_paths_3d_fast(int fr_total, 
	int fr_es, std::vector<double> parameters) {
	
	this->fr_total = fr_total; 
	this->fr_es = fr_es; 

	std::cout << "fr total: " << fr_total << ", fr es: " << fr_es << std::endl;

	// initialize frame weights 
	double a = .7;
	std::vector<double> fr_wts(fr_total, 1);
	for (int i = 0; i < fr_total; i++) {
		fr_wts[i] = 1/(1 + exp(-a*(i+1 - fr_es+2)));
	}
	this->fr_wts = fr_wts;

	// initialize optimization parameters by loading 
	DIST_IN_WT 		= parameters[0];
	DIST_PREV_WT 	= parameters[1];
	PT_DERIV_WT 	= parameters[2];
	SH_PREV_WT 		= parameters[3];

	std::cout << "DIST_IN_WT: "  << DIST_IN_WT << ", DIST_PREV_WT: " <<  DIST_PREV_WT 
	<< ", PT_DERIV_WT: " << PT_DERIV_WT << ", SH_PREV_WT: " << SH_PREV_WT  << std::endl; 
}



/**
 * @abstract - This method calculates the shortest path between
 * 		all given indices and end - normal
 * @param x - matrix  object with x points
 * @param y - matrix  object with y points
 * @param z - matrix  object with z points
 * @param near_nbor_inds - nearest neighbor indices for each point
 * 				each row has the neighbors, row accessed by (j*colSize + i)
 * @param paths - the return variable, the shortest path
 * @param path_dist - to store the path distance
 */ 
void shortest_paths_3d_fast::findShortestPaths(matrix<double> &x, matrix<double> y, matrix<double> &z, 
	matrix<double> &init_xyz, matrix3d &near_nbor_inds,
	matrix3d &nbors_dist_Mat, matrix3d &nbors_shp_Mat,
	std::vector<int> &inds_to_track, std::vector<int> &init_ref_inds, 
	matrix<int> &path_inds, std::vector<double> &path_dist) {
	
	std::cout << "using the 2nd order method" << std::endl;
    
	int num_paths = inds_to_track.size();

	double rand_large_val 	= 100000000;
	int n_states 			= x.size1(); // this is rowSize
	int n_pts 				= x.size2(); // this is colSize (points at each state)

	int nbors_length = near_nbor_inds.depth();
    
     // initialize required variables/matrices
    std::vector<double> large_pt(3, rand_large_val);
    std::vector<double> root_root_pt(3);

	std::vector<double> root_pt(3);
	std::vector<double> initial_pt(3);
	std::vector<double> current_pt(3);
	// std::vector<int> path_inds(n_states);

	// initializing nbor related vectors and such here
	std::vector<double> nnbors = std::vector<double>(nbors_length);
	std::vector<double> root_to_pt_euc_dist_vec = std::vector<double>(nbors_length);
	std::vector<double> root_to_pt_shp_dist_vec = std::vector<double>(nbors_length);
	double root_to_pt_euc_dist, root_to_pt_shp_dist;

	// initialize the dist_to and edge_to matrices
	matrix<double> dist_to = matrix<double>(n_states, n_pts);
	matrix<double> edge_to = matrix<double>(n_states, n_pts);
    
	// for each path to find
	for (int p_count = 0; p_count < num_paths; p_count++) {
		// assign large value to all points
		for (int i = 1; i < n_states; i++) {
			for (int j = 0; j < n_pts; j++) {
				dist_to(i, j) = rand_large_val;
			}
		}

		// row 1 initialize as zero
		for (int j = 0; j < n_pts; j++) {
			dist_to(0, j) = 0;
		}
		//-------------------------------------------//
		// std::cout << "before in sh " << std::endl;
		initial_pt[0] = init_xyz(0, init_ref_inds[p_count]);
		initial_pt[1] = init_xyz(1, init_ref_inds[p_count]);
		initial_pt[2] = init_xyz(2, init_ref_inds[p_count]);
		// std::cout << "after in sh " << std::endl;
       
		//-------------------------------------------//
		std::vector<int> roots_to_visit;
		roots_to_visit.push_back((int)inds_to_track[p_count]);

		// now update all remaining rows and columns
		for (int fr_no = 0; fr_no < n_states - 1; fr_no++) {
			set visited_branches = set();
			// std::cout << "visiting " << roots_to_visit.size() << " nodes in frame " << fr_no << std::endl; 
			for (int vi = 0; vi < roots_to_visit.size(); vi++) {
				int pt_no = roots_to_visit[vi];
				root_pt[0] = x(fr_no, pt_no);
				root_pt[1] = y(fr_no, pt_no);
				root_pt[2] = z(fr_no, pt_no);
				// std::cout << "before root sh " << std::endl;
                // std::cout << "before root sh " << std::endl;

                if (fr_no == 0) {
                	root_root_pt[0] = 0;
                	root_root_pt[1] = 0;
                	root_root_pt[2] = 0;
                } else {
                	// std::cout << "edge non zero: " << std::endl;
                	// std::cout << edge_to(fr_no, pt_no) << std::endl;	
                    root_root_pt[0] = x(fr_no - 1, edge_to(fr_no, pt_no));
                    root_root_pt[1] = y(fr_no - 1, edge_to(fr_no, pt_no));
                    root_root_pt[2] = z(fr_no - 1, edge_to(fr_no, pt_no));
                }
				// nbors
                // std::cout << "before nbors " << std::endl;                    
				nnbors = near_nbor_inds.at(fr_no, pt_no);	
				root_to_pt_euc_dist_vec = nbors_dist_Mat.at(fr_no, pt_no);
				root_to_pt_shp_dist_vec = nbors_shp_Mat.at(fr_no, pt_no);

                // std::cout << "after nbors " << std::endl;

				for (int j = 0; j < nbors_length; j++) {
					int nbor_j = (int) nnbors[j];
					// std::cout << "visiting neighbor ind: " << nbor_j << std::endl;
					current_pt[0] = x(fr_no + 1, nbor_j);
					current_pt[1] = y(fr_no + 1, nbor_j);
					current_pt[2] = z(fr_no + 1, nbor_j);

					root_to_pt_euc_dist = root_to_pt_euc_dist_vec[j];
					root_to_pt_shp_dist = root_to_pt_shp_dist_vec[j];


                    double dist_updated = dist_new(current_pt, root_pt, initial_pt, 
                    	root_root_pt, root_to_pt_euc_dist, root_to_pt_shp_dist, fr_no + 1);

					// if dist_new + dist_root less than dist old, update
					if (dist_updated + dist_to(fr_no, pt_no) < dist_to(fr_no + 1, nbor_j)) {
						/*std::cout << fr_no << ", " << pt_no << ", "<< "we updating: " <<  
						dist_to.at(fr_no + 1, pt_no) << ", to: " << dist_new + dist_to.at(fr_no, pt_no) << std::endl;*/
						dist_to(fr_no + 1, nbor_j) = dist_updated + dist_to(fr_no, pt_no);
						edge_to(fr_no + 1, nbor_j) = pt_no;
						visited_branches.insert(nbor_j);			
					}
				}
			}
			roots_to_visit = visited_branches.get_elements();		
			// std::cout << "how we doin' ?" <<  "to visit next frame: " << roots_to_visit.size() << std::endl;
		}

		//-------------------------------------------//
		// find path based on dist_to final row
		//std::cout << "min dist: " << info.val << ", ind: " << info.j << std::endl;
		int min_ind 		= minPosInRow(dist_to, n_states - 1);
        // min_ind             = inds_to_track[p_count];
        path_dist[p_count]  = dist_to(n_states - 1, min_ind);
		path_inds(p_count, n_states-1) = min_ind;

		// std::cout << "min dist: " << path_dist[p_count] << ", start ind: " << inds_to_track[p_count] << ", min ind: " << min_ind << "  ";

		for (int pt_no = n_states - 1; pt_no > 1; pt_no--) {
			path_inds(p_count, pt_no -1) = edge_to(pt_no, path_inds(p_count, pt_no));
		}
		// assign the 1st path points as the given path inds
		path_inds(p_count, 0) = inds_to_track[p_count];
		//-------------------------------------------//
        
        // update every 100th path no
        if (p_count % 100 == 0) {
	        std::cout << "path no: " << p_count << " done" << std::endl;        	
        }
	}
}



/**
 * find minimum, for a column, of a matrix
 * @param mat - matrix m
 * @param col - the column where minimum is being searched
 */
int shortest_paths_3d_fast::minPosInRow(matrix<double> &mat, int row) {
	double min_val = 1.7976931348623157e+308;
	int minJ = 0;
	for (int j = 0; j < mat.size2(); j++) {
		if (mat(row, j) < min_val) {
			min_val = mat(row, j);
			minJ 	= j;
		}
	}
	return minJ;
}


/**
 * multiply vector a by double b
 * @param a - double vector
 * @param b - double value
 */
std::vector<double> shortest_paths_3d_fast::multiply(std::vector<double> &a, double b) {
	std::vector<double> mult(a.size());
	for (int i = 0; i < a.size(); i++) {
		mult[i] = b*a[i];
	}
    return mult;
}


/**
 * subtract vector b from a
 * @param a - double vector 1
 * @param b - double vector 2
 */
std::vector<double> shortest_paths_3d_fast::subtract(std::vector<double> &a, std::vector<double> &b) {
	if (a.size() != b.size())
		throw std::invalid_argument("a, b sizes don't match");
	std::vector<double> diff(a.size());
	for (int i = 0; i < a.size(); i++) {
		diff[i] = a[i] - b[i];
	}
    return diff;
}


/**
 * add vector b and a
 * @param a - double vector 1
 * @param b - double vector 2
 */
std::vector<double> shortest_paths_3d_fast::add(std::vector<double> &a, std::vector<double> &b) {
	if (a.size() != b.size())
		throw std::invalid_argument("a, b sizes don't match");
	std::vector<double> sum(a.size());
	for (int i = 0; i < a.size(); i++) {
		sum[i] = a[i] + b[i];
	}
    return sum;
}


/**
 * distance between two vectors
 * @param a - double vector 1
 * @param b - double vector 2
 */
double shortest_paths_3d_fast::norm(std::vector<double> &a, std::vector<double> &b) {
	if (a.size() != b.size())
		throw std::invalid_argument("a, b sizes don't match");
	double dist = 0;
	for (int i = 0; i < a.size(); i++) {
		dist = dist + pow(a[i] - b[i], 2);
	}
    return sqrt(dist);
}




/**
 * distance metric according to different rule
 * @param current_pt - xyz location of node
 * @param root_pt - location of previous node
 * @param initial_pt - the xyz location of start node
 * @param current_shape - shape feature matrix at current node
 * @param root_shape - shape feature matrix at previous node
 * @param initial_shape - shape feature matrix at start node
 * @param root_root_pt - location of the root of root point
 * @param root_root_sh - shape of the root of the root point
 * @param fr_i - current frame of operation
 * @param fr_total - total frames
 */
double shortest_paths_3d_fast::dist_new(std::vector<double> &current_pt, std::vector<double> &root_pt,
			std::vector<double> &initial_pt, std::vector<double> &root_root_pt,
			double root_to_pt_euc_dist, double root_to_pt_shp_dist,  int fr_i) {

	/*
		DIST_IN_WT 		= parameters[0];
	DIST_PREV_WT 	= parameters[1];
	PT_DERIV_WT 	= parameters[2];
	SH_IN_WT 		= parameters[3];
	SH_PREV_WT 		= parameters[4];
	SH_DERIV_WT 	= parameters[5];
	*/

    // weights
    double FR_WT = fr_wts[fr_i];
    // double FR_WT = 1;

    double pt_dist_init 	= norm(current_pt, initial_pt);
    double pt_dist_prev 	= root_to_pt_euc_dist;
    double sh_dist_prev 	= root_to_pt_shp_dist;

    // std::cout << pt_dist_prev << std::endl;
    
    std::vector<double> sum_pt = add(current_pt, root_root_pt);
    
    std::vector<double> mult_pt = multiply(root_pt, 2);
    
    double pt_deriv = norm(sum_pt, mult_pt);
    

    double dist_val = (FR_WT*DIST_IN_WT*pt_dist_init + 
        DIST_PREV_WT*pt_dist_prev +
        PT_DERIV_WT*pt_deriv + 
        SH_PREV_WT*sh_dist_prev)/1000;

    return dist_val;
}




shortest_paths_3d_fast::~shortest_paths_3d_fast() {
	// TODO Auto-generated destructor stub
}

