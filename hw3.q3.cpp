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


void HW2_3(double xmin, double xmax, double ymin, double ymax, double tf){

    int N = 129;
    Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
    double dx = newGrid.get_dx();
    std::cout << "space step size : " << dx << std::endl;
    //    initial/exact solution:
    std::vector<double> ini_soln;
    ini_soln.assign(N*N, 0.);
    std::vector<double> level_set;
    level_set.assign(N*N, 0.);
    for (int n = 0; n < N * N; n++)
    {
        double x = newGrid.x_from_n(n);
        double y = newGrid.y_from_n(n);
        ini_soln[n] = ini_cond(x,y);
        level_set[n] = signum(ini_soln[n]);
    }
    char file_name[250];
    sprintf(file_name,"../q3VTK.vtk");
    newGrid.print_VTK_format(file_name);
    newGrid.print_VTK_format(ini_soln, "ini_soln", file_name);
    newGrid.print_VTK_format(level_set, "level_set_ini", file_name);

    SL_method SL(newGrid, ini_soln);
    std::vector<double> sol= ini_soln;
    double dt = dx;
    double t = 0.;
    int num_iter = 1;

    while(t < tf){
        dt = std::min(dt, tf-t);
        SL.set_initial(SL.get_sol());
        SL.advection_solver(dt);
        for (int i = 0; i < 10; i++){
            SL.godunov(dt/50.);
        }
        sol = SL.get_sol();
        char data_name[250];
        if (num_iter % 100 == 0 ){
            sprintf(data_name,"num_iter=%d", num_iter);
            newGrid.print_VTK_format(sol, data_name,file_name);
        }
        num_iter += 1;
        t += dt;
    }
    char data_name[250];
    sprintf(data_name,"final_soln");
    newGrid.print_VTK_format(sol, data_name,file_name);
    for (int n = 0; n < N * N; n++)
    {
        level_set[n] = signum(sol[n]);

    }

    newGrid.print_VTK_format(level_set, "final_level_set",file_name);



}