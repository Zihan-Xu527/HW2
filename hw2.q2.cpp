//
// Created by Zihan Xu on 10/23/22.
//

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "Grid2d.h"
#include "math_tools.h"
#include "SL_method.h"
#include "omp.h"


void HW2_2(double xmin, double xmax, double ymin, double ymax, double tf){

    int N = 129;
    Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
    double dx = newGrid.get_dx();
    std::cout << "space step size : " << dx << std::endl;

    //    initial solution:
    std::vector<double> ini_soln;
    ini_soln.assign(N*N, 0.);
    for (int n = 0; n < N * N; n++)
    {
        double x = newGrid.x_from_n(n);
        double y = newGrid.y_from_n(n);
        ini_soln[n] = ini_cond(x,y) * 0.5;
    }
    std::vector<double> true_soln;
    true_soln.assign(N*N, 0.);
    for (int n = 0; n < N * N; n++)
    {
        double x = newGrid.x_from_n(n);
        double y = newGrid.y_from_n(n);
        true_soln[n] = ini_cond(x,y) ;
    }
    newGrid.print_VTK_format("../q2reinit.vtk");
    newGrid.print_VTK_format(true_soln, "exact_soln","../q2reinit.vtk");

    //  REINITIALIZATION
    SL_method SL(newGrid, ini_soln);
    std::vector<double> sol= ini_soln;
    double d_tau = .5 * dx;
    double t = 0.;
//    int max_iter = 1000;
    int num_iter = 0;

    std::vector<double> diff;
    diff.assign(N*N, 0.);
    std::vector<double> err= err_norm(true_soln, ini_soln, 0.2, diff);
    double linf_err = err[2];


    while( t < tf &&  linf_err > 0.032){

        char data_name[250];
        if (num_iter % 100 == 0){
            sprintf(data_name,"num_iter=%d", num_iter);
            newGrid.print_VTK_format(sol, data_name,"../q2reinit.vtk");
            std::cout<<"max norm of iteration "<<num_iter<<": "<<linf_err<<std::endl;
        }
        d_tau = std::min(d_tau, tf-t);
        SL.godunov(d_tau);
        sol = SL.get_sol();
        err = err_norm(true_soln, sol, 0.1, diff);
        linf_err = err[2];
        num_iter += 1;

        t += d_tau;
    }
    std::vector<double> final_soln = SL.get_sol();
    err = err_norm(true_soln, sol, 0.1, diff);
    linf_err = err[2];
    std::cout <<"number of iterations: "<< num_iter << std::endl;
    std::cout <<"final err: "<< linf_err << std::endl;
    std::cout <<"final time: "<< t << std::endl;
    newGrid.print_VTK_format(final_soln, "final_soln","../q2reinit.vtk");
    newGrid.print_VTK_format(diff, "err","../q2reinit.vtk");
}