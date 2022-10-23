//
// Created by Zihan Xu on 10/22/22.
//

#include "SL_method.h"
#include <vector>

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