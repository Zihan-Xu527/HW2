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

    long N = 3;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    Grid2d newGrid(N, N, xmin, xmax, ymin, ymax);
    double dx = newGrid.get_dx();
    double dt = dx / 0.5; // various dx/dt ratios: 0.5, 1, 5, 10

    std::vector<double> ini_sol;
    ini_sol.resize(N*N);
    ini_sol.assign(ini_sol.size(), 0.);

    for (int i; i < N*N; i++)
    {
        double x = newGrid.x_from_n(i);
        double y = newGrid.y_from_n(i);
        ini_sol[i] = ini_cond(x,y);
    }

    ENO_interpolation(newGrid, ini_sol, -1., -1.);





//    std::cout << "This is dx: " << newGrid.get_dx() << std::endl;
//    std::cout << "This is dy: " << newGrid.get_dy() << std::endl;
//
//    std::cout << "Test i_from_n(7): " << newGrid.i_from_n(7) << std::endl;
//    std::cout << "Test j_from_n(7): " << newGrid.j_from_n(7) << std::endl;
//    std::cout << "Test n_from_ij(2,2): " << newGrid.n_from_ij(2, 2) << std::endl;
//
//    std::cout << "Test x_from_n(2): " << newGrid.x_from_n(2) << std::endl;
//    std::cout << "Test y_from_n(2): " << newGrid.y_from_n(2) << std::endl;
//
//    std::vector<double> values;
//    values.resize(N*M);
//    for (int i = 0; i < N*M; i++){
//        values[i] = i;
//    }
//
//    newGrid.print_VTK_format("../newVTK.vtk");
//    newGrid.print_VTK_format(values, "newData", "../newVTK.vtk");


    return 0;

}

