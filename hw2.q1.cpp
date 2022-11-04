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



void check_dx(){
    double tf = 2 * M_PI;
    std::ofstream output4;
    std::string fileName4 = "../check_dx.csv";
    output4.open(fileName4);
    // dx/dt ratios: 0.5, 1, 2, 4
    // dx : 1/64, 1/32, 1/16, 1/8
    // N : 129, 65, 33, 17
    int NN[4] = {17, 33, 65, 129};
    for (int i = 0; i < 4; i++){
        double dt = 1./16.; //fix time step size dt
        int N = NN[i] ;
        Grid2d newGrid(N, N, -1., 1., -1., 1.);
        double dx = newGrid.get_dx();
        output4 << dx << ",";
        std::vector<double> ini_soln;
        ini_soln.assign(N*N, 0.);
        for (int n = 0; n < (newGrid.get_N() * newGrid.get_M()) ; n++)
        {
            double x = newGrid.x_from_n(n);
            double y = newGrid.y_from_n(n);
            ini_soln[n] = ini_cond(x,y);
        }
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
        std::vector<double> err = err_norm(ini_soln, numer_soln, diff);
        output4 << err[0] << "," << err[1] << "," << err[2]<<"\n";
        std::cout<<"errors for dx="<<dx<<": "<<err[2]<<std::endl;
    }
    output4.close();

}



void HW2_1(double xmin, double xmax, double ymin, double ymax, double tf){

    int NN[4] = {17, 33, 65, 129};
//    int NN[4] = {9,17, 33, 65};
    double ratios[4] = {0.5, 1., 5., 10.}; // dx/dt ratios: 0.5, 1, 5, 10
    //output
    std::ofstream output1;
    std::string fileName1 = "../l1_errors.csv";
    output1.open(fileName1);

    std::ofstream output2;
    std::string fileName2 = "../l2_errors.csv";
    output2.open(fileName2);

    std::ofstream output3;
    std::string fileName3 = "../inf_errors.csv";
    output3.open(fileName3);

    output1<<" ,"<<ratios[0]<<","<<ratios[1]<<","<<ratios[2]<<","<<ratios[3]<<", "<<"\n";
    output2<<" ,"<<ratios[0]<<","<<ratios[1]<<","<<ratios[2]<<","<<ratios[3]<<", "<<"\n";
    output3<<" ,"<<ratios[0]<<","<<ratios[1]<<","<<ratios[2]<<","<<ratios[3]<<", "<<"\n";

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
            std::vector<double> err = err_norm(ini_soln, numer_soln, diff);
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
        // check convergence behavior for dt
        double dt = dx * dx;
        std::cout << "time step size dt = dx^2: " << dt << std::endl;

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
        std::vector<double> err = err_norm(ini_soln, numer_soln, diff);
        std::cout << "L1 norm: " << err[0] << std::endl;
        std::cout << "L2 norm: " << err[1] << std::endl;
        std::cout << "Max norm: " << err[2] << std::endl;
        output1<< err[0] ;
        output2<< err[1] ;
        output3<< err[2] ;

        output1 << "\n";
        output2 << "\n";
        output3 << "\n";

        std::cout <<"------------------------------------------------------------------------------------"<<std::endl;

    }
    output1.close();
    output2.close();
    output3.close();

//  if uncomment check_dx() function, we can use check_dx.csv to analyze the convergence behavior for dx
    check_dx();

}