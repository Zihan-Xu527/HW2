//
// Created by Zihan Xu on 10/22/22.
//

#ifndef F22_MATH_233_LABS_SL_METHOD_H
#define F22_MATH_233_LABS_SL_METHOD_H

#include "Grid2d.h"
#include <vector>
#include "math_tools.h"
#include <math.h>

// Semi-Langrangian Method
class SL_method {
private:
    Grid2d sl_grid;
    std::vector<double> ini_sol;
    std::vector<double> sol;
    double vel_x;
    double vel_y;
//    std::vector<double> vel_u;
//    std::vector<double> vel_v;
//    void find_trajectory(int n, double & x_d, double & y_d, double dt);

public:
    SL_method();
    SL_method(Grid2d grid, std::vector<double> ini_sol);
    void set_grid(Grid2d & new_grid){sl_grid = new_grid;};// set grid
    void set_init(std::vector<double> sol){ini_sol = sol;}; //set initial
    void set_velocity(double x, double y);
    std::vector<double> trajectory_interpolation(std::vector<double> & func, int n, double dt);
    void advection_solver(double dt);
    std::vector<double> get_sol(){ return sol; };        // access solution
    void reinitialize(double dt);
};


#endif //HW2_SL_METHOD_H
