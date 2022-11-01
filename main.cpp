#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "Grid2d.h"
#include "math_tools.h"
#include "SL_method.h"
#include "omp.h"



int main() {

    void HW2_1(double xmin, double xmax, double ymin, double ymax, double tf);
    void HW2_2(double xmin, double xmax, double ymin, double ymax, double tf);
    void HW2_3(double xmin, double xmax, double ymin, double ymax, double tf);
    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;
    double tf = 2*M_PI;

//    std::cout << "SOLUTION OF QUESTION 1: "<< std::endl;
//    HW2_1(xmin, xmax, ymin, ymax, tf);
//    std::cout << "SOLUTION OF QUESTION 2: "<< std::endl;
//    HW2_2(xmin, xmax, ymin, ymax, tf);
    std::cout << "SOLUTION OF QUESTION 3: "<< std::endl;
    HW2_3(xmin, xmax, ymin, ymax, tf);



    return 0;

}

