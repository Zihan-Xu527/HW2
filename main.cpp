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

int main() {
//    std::cout << "Hello, MATH233! :)" << std::endl;

    long N = 11;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;
    double tf = 2*M_PI;


    Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
    double dx = newGrid.get_dx();
    double dt = dx / 10.0; // various dx/dt ratios: 0.5, 1, 5, 10

    std::cout << "space step size : " << dx << std::endl;
    std::cout << "time step size : " << dt << std::endl;

    std::vector<double> ini_sol;
    ini_sol.resize(N*N);

//    initial solution:
    for (int i = 0; i < N*N; i++)
    {
        double x = newGrid.x_from_n(i);
        double y = newGrid.y_from_n(i);
        ini_sol[i] = ini_cond(x,y);

    }


    SL_method SL(newGrid, ini_sol);
    int num_iters = tf/dt;
    double dt_mod = tf - num_iters*dt;
    std::cout << "number of iterations: " << num_iters << std::endl;

    std::vector<double> true_soln;
    true_soln.resize(N*N);
    for(int n=0 ; n < (N*N) ; n++) {
        double x = newGrid.x_from_n(n);
        double y = newGrid.y_from_n(n);
        true_soln[n] = std::sqrt( pow((x * cos(tf) + y * sin(tf) - 0.25) , 2)
                             +  pow ((y * cos(tf) - x * sin(tf)) ,2)) - 0.2;

    }

    std::vector<double> numer_soln;
    numer_soln.resize(N*N);
    for (int i=0; i < num_iters; i++){
        SL.set_init(SL.get_sol());
        SL.advection_solver(dt);
    }
    SL.set_init(SL.get_sol());
    SL.advection_solver(dt_mod);
    numer_soln = SL.get_sol();

    std::vector<double> err;
    err.resize(N*N);
    double sum = 0.;
    double sum2 = 0.;

    for (int i = 0; i<N*N; i++){
        err[i] = abs(numer_soln[i]-true_soln[i]);
        sum += err[i];
        sum2 += err[i]*err[i];
    }
    double max_err = *max_element(err.begin(), err.end());
    std::cout << "L1 norm: " << sum << std::endl;
    std::cout << "L2 norm: " << std::sqrt(sum2) << std::endl;
    std::cout << "Max norm: "<< max_err << std::endl;




    return 0;

}

