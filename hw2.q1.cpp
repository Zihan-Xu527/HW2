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
#include <fstream>


std::vector<double> err_norm(std::vector<double> x, std::vector<double> y, double eps, std::vector<double> & diff){
    std::vector<double> err;
    err.resize(3); // reserve l1, l2, and max norm respectively
    double max = 0.;
//#pragma omp parallel for
    for (int i = 0; i < x.size(); i++){
        diff[i] = std::abs(x[i] - y[i]);
        err[0] += diff[i]; //l1 norm
        err[1] += pow(diff[i], 2);
        max = diff[i] > max ? diff[i] : max; //max norm
//        if (std::abs(x[i]) < eps){
//            err[0] += diff[i]; //l1 norm
//            err[1] += pow(diff[i], 2);
//            max = diff[i] > max ? diff[i] : max; //max norm
//        }

    }
    err[1] = std::sqrt(err[1]); //l2 norm
    err[2] = max;
    return err;
}

void HW2_1(double xmin, double xmax, double ymin, double ymax, double tf){

    int NN[4] = {17, 33, 65, 129};
    double ratios[4] = {0.5, 1., 5., 10.}; // dx/dt ratios: 0.5, 1, 5, 10
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

    for (int i = 0; i < 4 ; i++){

        int N = NN[i];
        std::cout << "N = " << N <<" , ";
        Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
        double dx = newGrid.get_dx();
        std::cout << "space step size : " << dx << std::endl;
        output1<< dx <<",";
        output2<< dx <<",";
        output3<< dx <<",";

        //    initial/exact solution:
        std::vector<double> ini_soln;
        ini_soln.assign(N*N, 0.);
        for (int n = 0; n < N * N; n++)
        {
            double x = newGrid.x_from_n(n);
            double y = newGrid.y_from_n(n);
            ini_soln[n] = ini_cond(x,y);
        }
        char file_name[250];
        sprintf(file_name,"../q1VTK%d.vtk", NN[i]);
        newGrid.print_VTK_format(file_name);
        newGrid.print_VTK_format(ini_soln, "iniData", file_name);


        for (int j = 0; j < 4 ; j++) {
            double ratio = ratios[j];
            std::cout << "dx/dt = " << ratio << ", ";
            double dt = dx/ratio;
            std::cout << "time step size : " << dt << std::endl;


            SL_method SL(newGrid, ini_soln);
            double t = 0.;
            std::vector<double> numer_soln;
            numer_soln.assign(ini_soln.size(), 0.);
            while(t < tf){
                dt = std::min(dt, tf-t);
                SL.advection_solver(dt);
                SL.set_initial(SL.get_sol());
                t += dt;
            }
            numer_soln = SL.get_sol();
            std::vector<double> diff;
            diff.assign(ini_soln.size(), 0.);
            std::vector<double> err = err_norm(ini_soln, numer_soln, 0.1, diff);
            std::cout << "L1 norm: " << err[0] << std::endl;
            std::cout << "L2 norm: " << err[1] << std::endl;
            std::cout << "Max norm: " << err[2] << std::endl;

            output1<< err[0] <<",";
            output2<< err[1] <<",";
            output3<< err[2] <<",";




            char numer_name[250];
            sprintf(numer_name, "numer_sol=%d", j);
            newGrid.print_VTK_format(numer_soln, numer_name, file_name);
            char err_name[250];
            sprintf(err_name, "errors=%d", j);
            newGrid.print_VTK_format(diff, err_name, file_name);

        }
        output1 << "\n";
        output2 << "\n";
        output3 << "\n";

        std::cout <<"------------------------------------------------------------------------------------"<<std::endl;

    }
    output1.close();
    output2.close();
    output3.close();

}