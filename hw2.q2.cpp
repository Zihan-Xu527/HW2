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

    //    initial/exact solution:
    std::vector<double> ini_soln;
    ini_soln.assign(N*N, 0.);
    for (int n = 0; n < N * N; n++)
    {
        double x = newGrid.x_from_n(n);
        double y = newGrid.y_from_n(n);
        ini_soln[n] = ini_cond(x,y)*0.5;
    }
    newGrid.print_VTK_format("../q2reinit.vtk");
    newGrid.print_VTK_format(ini_soln, "ini_soln","../q2reinit.vtk");

    //  REINITIALIZATION
    SL_method SL(newGrid, ini_soln);
    std::vector<double> sol= ini_soln;
//    sol.resize(N*N);
    double d_tau = .5 * dx;
    double t = 0.;
    int num_iter = 0;

//    int num_iter = floor(tf/d_tau);

//    for(int i = 0; i < num_iter; i++){
//        char data_name[250];
//        if (i%100 == 0){
//            sprintf(data_name,"num_iter=%d", i);
//            newGrid.print_VTK_format(ini_soln, data_name,"reinit.vtk");
//        }
//        SL.godunov(d_tau);
//        sol = SL.get_sol();
//
//
//    }


    while( t < tf ){
        num_iter += 1;
        d_tau = std::min(d_tau, tf-t);
        SL.godunov(d_tau);
        sol = SL.get_sol();
        char data_name[250];
        if (num_iter%100 == 0){
            sprintf(data_name,"num_iter=%d", num_iter);
            newGrid.print_VTK_format(sol, data_name,"../q2reinit.vtk");
        }
        t += d_tau;
    }
    std::vector<double> final_soln = SL.get_sol();
//    std::cout <<"dt = " <<d_tau<<std::endl;
    std::cout <<"number of iterations: "<< num_iter << std::endl;
    std::cout <<"final time: "<< t << std::endl;
    newGrid.print_VTK_format(final_soln, "final_soln","../q2reinit.vtk");


}