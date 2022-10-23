//
// Created by Zihan Xu on 10/20/22.
//

#ifndef HW2_MATH_TOOLS_H
#define HW2_MATH_TOOLS_H
#include "Grid2d.h"
#include <vector>

double bilinear_interpolation(Grid2d & grid, std::vector<double> & func, const double x, const double y);
double ENO_interpolation(Grid2d & grid, std::vector<double> & func, const double x, const double y);
int MAX(int a, int b);
double MAX(double a, double b);
double minmod(double a, double b);
double central_diff(double lo, double mid, double hi, double dx);
double sec_der_dx(Grid2d & grid, std::vector<double> & func, int n);
double sec_der_dy(Grid2d & grid, std::vector<double> & func, int n);
double bwd_dx(Grid2d & grid, std::vector<double> & func, int n);
double fwd_dx(Grid2d & grid, std::vector<double> & func, int n);
double bwd_dy(Grid2d & grid, std::vector<double> & func, int n);
double fwd_dy(Grid2d & grid, std::vector<double> & func, int n);
double signfunc(double x);
#endif //HW2_MATH_TOOLS_H
