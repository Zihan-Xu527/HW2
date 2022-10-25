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
//#include "omp.h"

double ini_cond(double x, double y){
    return std::sqrt( pow((x-0.25), 2) + y*y ) - 0.2;
}

std::vector<double> err_norm(std::vector<double> & x, std::vector<double> & y, double eps, std::vector<double> & diff){
    std::vector<double> err;
    err.resize(3); // reserve l1, l2, and max norm respectively

    for (int i = 0; i < x.size(); i++){
        diff[i] = std::abs(x[i] - y[i]);
        err[1] += pow(diff[i], 2);
        if (std::abs(x[i]) < eps)
            err[0] += diff[i]; //l1 norm
    }
    err[1] = std::sqrt(err[1]); //l2 norm
    err[2] = *max_element(diff.begin(), diff.end()); //max norm

    return err;
}
void HW2_1(double xmin, double xmax, double ymin, double ymax, double tf){




    int NN[4] = {9, 17, 33, 65};
    for (int i = 0; i < 4; i++){
        int N = NN[i];
        std::cout << "N = " << N <<" , ";
        Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
        double dx = newGrid.get_dx();
        std::cout << "space step size : " << dx << std::endl;

        //    initial/exact solution:
        std::vector<double> ini_sol;
        ini_sol.resize(N*N);
        for (int i = 0; i < N*N; i++)
        {
            double x = newGrid.x_from_n(i);
            double y = newGrid.y_from_n(i);
            ini_sol[i] = ini_cond(x,y);

        }
        char file_name[250];
        sprintf(file_name,"../allVTK%d.vtk", NN[i]);
        newGrid.print_VTK_format(file_name);
        newGrid.print_VTK_format(ini_sol, "iniData", file_name);


        double ratios[4] = {0.5, 1., 5., 10.};
        for (int j = 0; j < 4; j++) {
            double ratio = ratios[j];
            std::cout << "dx/dt = " << ratio << ", ";
            double dt = dx/ratio;
            std::cout << "time step size : " << dt << std::endl;
            int num_iters = floor(tf/dt);
            double dt_mod = tf - num_iters * dt;
            std::cout << "number of iterations: " << num_iters + 1 << std::endl;

            SL_method SL(newGrid, ini_sol);
            std::vector<double> numer_sol;
            numer_sol.resize(N*N);
            for (int k = 0; k < num_iters; k++) {
                numer_sol = SL.get_sol();
                SL.set_init(numer_sol);
                SL.advection_solver(dt);
//                SL.set_init(SL.get_sol());
            }
            SL.set_init(SL.get_sol());
            SL.advection_solver(dt_mod);
            numer_sol = SL.get_sol();

            std::vector<double> diff;
            diff.assign(ini_sol.size(), 0.);
            std::vector<double> err = err_norm(ini_sol, numer_sol, 0.1, diff);
            std::cout << "L1 norm: " << err[0] << std::endl;
            std::cout << "L2 norm: " << err[1] << std::endl;
            std::cout << "Max norm: " << err[2] << std::endl;




            char numer_name[250];
            sprintf(numer_name, "numer_sol=%d", j);
            newGrid.print_VTK_format(numer_sol, numer_name, file_name);
            char err_name[250];
            sprintf(err_name, "errors=%d", j);
            newGrid.print_VTK_format(diff, err_name, file_name);

        }

    }


}


