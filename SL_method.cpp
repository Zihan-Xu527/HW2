//
// Created by Zihan Xu on 10/22/22.
//

#include "SL_method.h"
#include <vector>
#include "math_tools.h"

SL_method::SL_method()
{

}
SL_method::SL_method(Grid2d grid, std::vector<double> ini_sol, double vel_x0, double vel_y0){
    sl_grid = grid;
    sol.resize(grid.get_N()*grid.get_M());
    sol = ini_sol;
    set_velocity(vel_x0,vel_y0);
    ini_sol.resize(grid.get_N()*grid.get_M());
    ini_sol = ini_sol;
}

void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
    // RK2 Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    set_velocity(x_0,y_0);
    double x_star = x_0 - dt/2 * vel_x;
    double y_star = y_0 - dt/2 * vel_y;
    set_velocity(x_star,y_star);
    x_d = x_0 - dt * vel_x;
    y_d = y_0 - dt * vel_y;
}

void SL_method::set_velocity(double vel_x0, double vel_y0){
    vel_x =  -vel_x0;;
    vel_y = vel_y0;
}

std::vector<double> SL_method::trajectory_interpolation(std::vector<double> &func, int n, double dt) {
    int N = sl_grid.get_N();
    int M = sl_grid.get_M();
    if ( int(func.size()) != (N*M) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");
    std::vector<double> coor;
    coor.resize(2);
    double x_d, y_d;

    find_trajectory(n, x_d, y_d, dt);
    coor[0] = x_d;
    coor[1] = y_d;
    return coor;
}