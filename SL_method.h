//
// Created by Zihan Xu on 10/22/22.
//

#ifndef F22_MATH_233_LABS_SL_METHOD_H
#define F22_MATH_233_LABS_SL_METHOD_H

#include "Grid2d.h"
#include <vector>

// Semi-Langrangian Method
class SL_method {
private:
    Grid2d sl_grid;
    std::vector<double> sol;
    double vel_x;
    double vel_y;
//    std::vector<double> vel_u;
//    std::vector<double> vel_v;
    void find_trajectory(int n, double & x_d, double & y_d, double dt);

public:
    void set_grid(Grid2d & new_grid){sl_grid = new_grid;} // set grid
    std::vector<double> get_sol(){ return sol; }        // access solution
    void set_velocity(double vel_u0, double vel_v0);
    //    void set_velocity(std::vector<double> & vel_u0, std::vector<double> & vel_v0);
//    std::vector<double> trajectory_interpolation(std::vector<double> & func, int n, double dt);


};


#endif //HW2_SL_METHOD_H
