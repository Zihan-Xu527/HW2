//
// Created by Zihan Xu on 10/19/22.
//

#ifndef LAB01_GRID2D_H
#define LAB01_GRID2D_H

#include <iostream>
#include <fstream>
#include <vector>

class Grid2d {
    // Objective: creates a 2D grid
private:
    double dx;  // spacing in x
    double dy;  // spacing in y
    long N;     // number of nodes in x
    long M;     // number of nodes in y
    double xmin;
    double xmax;
    double ymin;
    double ymax;
public:
    double get_dx();    // return dx to user
    double get_dy();    // return dy to user


    double get_xmin();
    double get_ymin();
    double get_xmax();
    double get_ymax();

    long get_N();
    long get_M();

    // Constructor to create object
    Grid2d();
    // Constructor to initialize values
    Grid2d(long NN, long MM, double xlo, double xhi, double ylo, double yhi);
    // Return coord (i,j) from grid index n
    long i_from_n(long n);
    long j_from_n(long n);
    // Return grid index n from coord (i, j)
    long n_from_ij(long i, long j);
    // Return position of grid index n
    double x_from_n(long n);
    double y_from_n(long n);

    // output file in VTK format
    void print_VTK_format(std::string output_file);
    void print_VTK_format(std::vector<double> &F, std::string data_name,
                          std::string file_name);

//    double fwd_dx(std::vector<double> &func, int n);
//    double bwd_dx(std::vector<double> &func, int n);
//
//    double fwd_dy(std::vector<double> &func, int n);
//    double bwd_dy(std::vector<double> &func, int n);
//
//    double sec_der_dx(std::vector<double> &func, int n);
//    double sec_der_dy(std::vector<double> &func, int n);
};


#endif //LAB01_GRID2D_H
