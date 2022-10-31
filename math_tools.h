//
// Created by Zihan Xu on 10/20/22.
//

#ifndef HW2_MATH_TOOLS_H
#define HW2_MATH_TOOLS_H
#include "Grid2d.h"
#include <vector>

double bilinear_interpolation(Grid2d & grid, std::vector<double> & func, double x, double y);
double minmod(double a, double b);
double central_diff(double lo, double mid, double hi, double dx);
double sec_der_dx(Grid2d & grid, std::vector<double> & func, int n);
double sec_der_dy(Grid2d & grid, std::vector<double> & func, int n);
double ENO_interpolation(Grid2d & grid, std::vector<double> & func, double x, double y);

double signum(double x);
double bwd_dx(Grid2d & grid, std::vector<double> & func, int n);
double fwd_dx(Grid2d & grid, std::vector<double> & func, int n);
double bwd_dy(Grid2d & grid, std::vector<double> & func, int n);
double fwd_dy(Grid2d & grid, std::vector<double> & func, int n);

double ini_cond(double x, double y);
#endif //HW2_MATH_TOOLS_H
