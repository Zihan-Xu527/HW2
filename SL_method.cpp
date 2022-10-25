//
// Created by Zihan Xu on 10/22/22.
//

#include "SL_method.h"
#include <vector>
#include "math_tools.h"
#include <math.h>

SL_method::SL_method()
{

}
SL_method::SL_method(Grid2d grid, std::vector<double> ini){
    sl_grid = grid;
    ini_sol.resize(grid.get_N()*grid.get_M());
    ini_sol = ini;
    sol.resize(grid.get_N()*grid.get_M());
    sol = ini_sol;
    vel_x = 0.;
    vel_y = 0.;

}

//void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
//    // RK2 Method
//
//
//    double x_0 = sl_grid.x_from_n(n);
//    double y_0 = sl_grid.y_from_n(n);
//    set_velocity(x_0,y_0);
//    double x_star = x_0 - dt/2 * vel_x;
//    double y_star = y_0 - dt/2 * vel_y;
//    set_velocity(x_star,y_star);
//    x_d = x_0 - dt * vel_x;
//    y_d = y_0 - dt * vel_y;
//}


void SL_method::set_velocity(double x, double y){
    vel_x =  -1.*y;
//    std::cout << "velocity of x: " << vel_x << std::endl;
    vel_y = x;
//    std::cout << "velocity of y: " << vel_y << std::endl;
}

//double SL_method::check_boundary(double x, double min, double max)
//{
//    if ( x < min )
//        return min;
//    else if ( x > max )
//        return max;
//    else
//        return x;
//
//}

std::vector<double> SL_method::trajectory_interpolation(std::vector<double> &func, int n, double dt) {
    int N = sl_grid.get_N();
    int M = sl_grid.get_M();
    if ( int(func.size()) != (N*M) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

        // RK2 Method


    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
//    std::cout << "x0: " << x_0 << ", y0: "<<y_0<<std::endl;

    set_velocity(x_0,y_0);
//    std::cout << "velocity of x: " << vel_x << ", velocity of y: "<<vel_y<<std::endl;

    double x_mid = x_0 - .5 * dt * vel_x;
    double y_mid = y_0 - .5 * dt * vel_y;
//    std::cout << "x_mid: " << x_mid << ", y_mid: "<<y_mid<<std::endl;
//    if (x_mid <= sl_grid.get_xmin()) { x_mid = sl_grid.get_xmin(); }
//    if (x_mid >= sl_grid.get_xmax()) { x_mid = sl_grid.get_xmax(); }
//
//    if (y_mid <= sl_grid.get_ymin()) { y_mid = sl_grid.get_ymin(); }
//    if (y_mid >= sl_grid.get_ymax()) { y_mid = sl_grid.get_ymax(); }

    set_velocity(x_mid,y_mid);
    double x_d = x_0 - dt * vel_x;
    double y_d = y_0 - dt * vel_y;
//    std::cout << "velocity of x: " << vel_x << ", velocity of y: "<<vel_y<<std::endl;
//    std::cout << "x_d: " << x_d << ", y_d: "<<y_d<<std::endl;


    std::vector<double> depart_coord;
    depart_coord.resize(2);

//    find_trajectory(n, x_d, y_d, dt);
    depart_coord[0] = x_d;
    depart_coord[1] = y_d;
    return depart_coord;
}

void SL_method::advection_solver(double dt)
{
    #pragma omp parallel for
    for (int n=0; n < (sl_grid.get_N()*sl_grid.get_M()) ; n++) {
        std::vector<double> depart_coord = trajectory_interpolation( ini_sol, n, dt );
        sol[n] = ENO_interpolation(sl_grid, ini_sol, depart_coord[0], depart_coord[1]);
    }

}

void SL_method::godunov(double dt){
    std::vector<double> func = sol;
    for (int n=0; n < (sl_grid.get_N()*sl_grid.get_M()) ; n++) {
        double phi_x = 0.;
        double phi_y = 0.;

        if (signfunc(ini_sol[n])* bwd_dx(sl_grid, func, n) <= 0. && signfunc(ini_sol[n])* fwd_dx(sl_grid, func, n) <= 0.)
            phi_x = fwd_dx(sl_grid, func, n);
        else if (signfunc(ini_sol[n])* bwd_dx(sl_grid, func, n) >= 0. && signfunc(ini_sol[n])* fwd_dx(sl_grid, func, n) >= 0.)
            phi_x = bwd_dx(sl_grid, func, n);
        else if (signfunc(ini_sol[n])* bwd_dx(sl_grid, func, n) >= 0. && signfunc(ini_sol[n])* fwd_dx(sl_grid, func, n) <= 0.)
            if ( abs(bwd_dx(sl_grid, func, n))>= abs(fwd_dx(sl_grid, func, n)))
                phi_x = bwd_dx(sl_grid, func, n);
            else
                phi_x = fwd_dx(sl_grid, func, n);

        if (signfunc(ini_sol[n])* bwd_dy(sl_grid, func, n) <= 0. && signfunc(ini_sol[n])* fwd_dy(sl_grid, func, n) <= 0.)
            phi_y = fwd_dy(sl_grid, func, n);
        else if (signfunc(ini_sol[n])* bwd_dy(sl_grid, func, n) >= 0. && signfunc(ini_sol[n])* fwd_dy(sl_grid, func, n) >= 0.)
            phi_y = bwd_dy(sl_grid, func, n);
        else if (signfunc(ini_sol[n])* bwd_dy(sl_grid, func, n) >= 0. && signfunc(ini_sol[n])* fwd_dy(sl_grid, func, n) <= 0.)
            if ( abs(bwd_dy(sl_grid, func, n))>= abs(fwd_dy(sl_grid, func, n)))
                phi_y = bwd_dy(sl_grid, func, n);
            else
                phi_y = fwd_dy(sl_grid, func, n);

        sol[n] = func[n] - dt * signfunc(ini_sol[n]) * ( sqrt(pow(phi_x,2)+pow(phi_y,2)) - 1. );
    }

}

