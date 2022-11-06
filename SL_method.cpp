//
// Created by Zihan Xu on 10/22/22.
//

#include "SL_method.h"
#include <vector>
#include "math_tools.h"
#include <math.h>
#include <cmath>

SL_method::SL_method(){

}

SL_method::SL_method(Grid2d grid, std::vector<double> ini_soln){
    sl_grid = grid;
    ini_sol.resize(grid.get_N()*grid.get_M());
    ini_sol = ini_soln;
    sol.resize(grid.get_N()*grid.get_M());
    sol = ini_sol;
    vel_x = 0.;
    vel_y = 0.;
}

void SL_method::set_velocity(const double x, const double y){
    vel_x = - 1. * y;
    vel_y = x;
}


std::vector<double> SL_method::trajectory_interpolation(std::vector<double> &func, int n, double dt) {

    if ( int(func.size()) != (sl_grid.get_N()*sl_grid.get_M()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

    // RK2 Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);

    set_velocity(x_0,y_0);
    double x_mid = x_0 - .5 * dt * vel_x;
    double y_mid = y_0 - .5 * dt * vel_y;

    set_velocity(x_mid,y_mid);
    double x_d = x_0 - dt * vel_x;
    double y_d = y_0 - dt * vel_y;

    std::vector<double> depart_coord{x_d, y_d};

    return depart_coord;
}

void SL_method::advection_solver(double dt)
{
//#pragma omp parallel for
    for (int n = 0 ; n < ( sl_grid.get_N() * sl_grid.get_M() ) ; n++) {
        std::vector<double> depart_coord = trajectory_interpolation( ini_sol, n, dt );
        sol[n] = ENO_interpolation(sl_grid, ini_sol, depart_coord[0], depart_coord[1]);
    }
}

void SL_method::godunov(double dt){
    std::vector<double> func = sol;
//#pragma omp parallel for
    for (int n=0; n < ( sl_grid.get_N() * sl_grid.get_M() ) ; n++) {
        double phi_x = 0.; // default condition c
        // condition a:
        if ( signum(ini_sol[n])* bwd_dx(sl_grid, func, n) <= 0. && signum(ini_sol[n])* fwd_dx(sl_grid, func, n) <= 0. )
        {
            phi_x = fwd_dx(sl_grid, func, n);
        }
        // condition b:
        if ( signum(ini_sol[n])* bwd_dx(sl_grid, func, n) >= 0. && signum(ini_sol[n])* fwd_dx(sl_grid, func, n) >= 0. )
        {
            phi_x = bwd_dx(sl_grid, func, n);
        }
        // condition d
        if ( signum(ini_sol[n])* bwd_dx(sl_grid, func, n) >= 0. && signum(ini_sol[n])* fwd_dx(sl_grid, func, n) <= 0. ){
            // d. i. or ii.
            phi_x = ( std::abs(bwd_dx(sl_grid, func, n))>= std::abs(fwd_dx(sl_grid, func, n)) ) ? bwd_dx(sl_grid, func, n) : fwd_dx(sl_grid, func, n);

        }



        double phi_y = 0.; // default condition c
        // condition a:
        if ( signum(ini_sol[n])* bwd_dy(sl_grid, func, n) <= 0. && signum(ini_sol[n])* fwd_dy(sl_grid, func, n) <= 0. )
        {
            phi_y = fwd_dy(sl_grid, func, n);
        }
        // condition b:
        if (signum(ini_sol[n])* bwd_dy(sl_grid, func, n) >= 0. && signum(ini_sol[n])* fwd_dy(sl_grid, func, n) >= 0.)
        {
            phi_y = bwd_dy(sl_grid, func, n);
        }
        // condition d
        if (signum(ini_sol[n])* bwd_dy(sl_grid, func, n) >= 0. && signum(ini_sol[n])* fwd_dy(sl_grid, func, n) <= 0.)
        {   // d. i. or ii.
            phi_y = ( std::abs(bwd_dy(sl_grid, func, n))>= std::abs(fwd_dy(sl_grid, func, n)) ) ? bwd_dy(sl_grid, func, n) : fwd_dy(sl_grid, func, n);
        }

        sol[n] = func[n] - dt * signum(ini_sol[n]) * ( std::sqrt( pow(phi_x,2) + pow(phi_y,2) ) - 1. );
    }

}

