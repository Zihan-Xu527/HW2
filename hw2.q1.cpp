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
#include <fstream>

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
    double ratios[4] = {0.5, 1., 5., 10.};
    //output
    std::ofstream output1;
    std::string fileName1 = "l1_errors.csv";
    output1.open(fileName1);

    std::ofstream output2;
    std::string fileName2 = "l2_errors.csv";
    output2.open(fileName2);

    std::ofstream output3;
    std::string fileName3 = "inf_errors.csv";
    output3.open(fileName3);




    output1<<" ,"<<ratios[0]<<","<<ratios[1]<<","<<ratios[2]<<","<<ratios[3]<<"\n";
    output2<<" ,"<<ratios[0]<<","<<ratios[1]<<","<<ratios[2]<<","<<ratios[3]<<"\n";
    output3<<" ,"<<ratios[0]<<","<<ratios[1]<<","<<ratios[2]<<","<<ratios[3]<<"\n";
    for (int i = 0; i < 4; i++){
        int N = NN[i];
        std::cout << "N = " << N <<" , ";
        Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
        double dx = newGrid.get_dx();
        std::cout << "space step size : " << dx << std::endl;
        output1<< dx <<",";
        output2<< dx <<",";
        output3<< dx <<",";

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

            output1<< err[0] <<",";
            output2<< err[1] <<",";
            output3<< err[2] <<",";




            char numer_name[250];
            sprintf(numer_name, "numer_sol=%d", j);
            newGrid.print_VTK_format(numer_sol, numer_name, file_name);
            char err_name[250];
            sprintf(err_name, "errors=%d", j);
            newGrid.print_VTK_format(diff, err_name, file_name);

        }
        output1 << "\n";
        output2 << "\n";
        output3 << "\n";

    }
    output1.close();
    output2.close();
    output3.close();


}


