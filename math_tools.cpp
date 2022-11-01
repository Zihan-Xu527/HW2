//
// Created by Zihan Xu on 10/20/22.
//

#include "math_tools.h"
#include <cmath>

// suppose we have grid in 2D [xmin,xmax] x [ymin,ymax]
// find cell in which (x,y) belongs
// find weighted avg of values (?)

double bilinear_interpolation(Grid2d & grid, std::vector<double> & func, double x, double y){

    double phi;
    double dx = grid.get_dx();
    double dy = grid.get_dy();
    double xmin = grid.get_xmin();
    double xmax = grid.get_xmax();
    double ymin = grid.get_ymin();
    double ymax = grid.get_ymax();
    int N = grid.get_N();
    int M = grid.get_M();
    int i,j;

    if ( int(func.size()) != (N*M) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

    // if (x,y) are outside domain, we need to find the nearest grid
    if (x >  xmax)
        i = N - 2;
    else if (x < xmin)
        i = 0;
    else
        i = floor((x - xmin)/dx);

    if (y > ymax)
        j = M - 2;
    else if (y < ymin)
        j = 0;
    else
        j = floor((y - ymin)/dy);

    double x_i = xmin + i * dx;
    double y_j = ymin + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    // Use quadratic interpolation (formula in Lab 3) to get value at x
    // (i.e. think weighted avg)
    // (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1) are the corners of the cell C
    phi  = func[grid.n_from_ij(i    ,   j  )]  * ( x_ip1 - x   ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1,   j  )]  * ( x     - x_i ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i    , j+1)]  * ( x_ip1 - x   ) * ( y     - y_j ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1, j+1)]  * ( x     - x_i ) * ( y     - y_j ) / (dx*dy) ;

    return phi;
}

double minmod(double a, double b){
    if ( a*b < 0.0 )
        return 0.;
    else if ( std::abs(a) < std::abs(b) )
        return a;
    else
        return b;
}

double central_diff(double lo, double mid, double hi, double h){
    return (lo - 2. * mid + hi) / (h * h);
}

double sec_der_dx(Grid2d & grid, std::vector<double> & func, int n)
{

    double dx = grid.get_dx();

    if ( grid.x_from_n(n) == grid.get_xmin() )
    {
        return central_diff(func[n], func[n + 1],func[n + 2], dx);
    }
    else if ( grid.x_from_n(n) == grid.get_xmax() )
    {
        return central_diff(func[n], func[n - 1],func[n - 2], dx);

    }
    else
    {
        return central_diff(func[n - 1], func[n],func[n + 1], dx);
    }
}

double sec_der_dy(Grid2d & grid, std::vector<double> & func, int n)
{
    int N = grid.get_N();
    double dy = grid.get_dy();

    if ( grid.y_from_n(n) == grid.get_ymin() )
    {
        return central_diff(func[n], func[n + N],func[n + 2 * N], dy);
    }
    else if ( grid.y_from_n(n) == grid.get_ymax() )
    {
        return central_diff(func[n], func[n - N],func[n - 2 * N], dy);

    }
    else
    {
        return central_diff(func[n - N], func[n],func[n + N], dy);
    }
}
// 01 --- 11
// |       |
// |       |
// 00 --- 10
double ENO_interpolation(Grid2d & grid, std::vector<double> & func, double x, double y){
    double phi;
    double dx = grid.get_dx();
    double dy = grid.get_dy();
    double xmin = grid.get_xmin();
    double xmax = grid.get_xmax();
    double ymin = grid.get_ymin();
    double ymax = grid.get_ymax();
    int N = grid.get_N();
    int M = grid.get_M();
    int i,j;

    if ( int(func.size()) != (N*M) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

    // if (x, y) is inside the domain, we take the floor of (x - xmin)/dx and (y - ymin)/dy
    // else if (x,y) is outside domain, we need to find the nearest grid

    if (x >= xmin && x <= xmax)
        i = floor((x - xmin)/dx);
    else
        i = x > xmax ? N - 2: 0;

    if (y >= ymin && y <= ymax)
        j = floor((y - ymin)/dy);
    else
        j = y > ymax ? M - 2: 0;

    double x_i = xmin + i * dx;
    double y_j = ymin + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    // Use quadratic interpolation (formula in Lab 3) to get value at x
    // (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1) are the corners of the cell C
    phi  = func[grid.n_from_ij(i    ,   j  )]  * ( x_ip1 - x   ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1,   j  )]  * ( x     - x_i ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i    , j+1)]  * ( x_ip1 - x   ) * ( y     - y_j ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1, j+1)]  * ( x     - x_i ) * ( y     - y_j ) / (dx*dy) ;

    double phi_xx_00 = sec_der_dx(grid, func, grid.n_from_ij(i    ,   j  ));
    double phi_xx_10 = sec_der_dx(grid, func, grid.n_from_ij(i+1,   j  ));
    double phi_xx_01 = sec_der_dx(grid, func, grid.n_from_ij(i    , j+1));
    double phi_xx_11 = sec_der_dx(grid, func, grid.n_from_ij(i+1, j+1));
    double phi_xx = minmod(minmod(phi_xx_00, phi_xx_01),minmod(phi_xx_10, phi_xx_11));
    phi -= .5 * (x - x_i) * (x_ip1 - x) * phi_xx;

    double phi_yy_00 = sec_der_dy(grid, func, grid.n_from_ij(i    ,   j  ));
    double phi_yy_10 = sec_der_dy(grid, func, grid.n_from_ij(i+1,   j  ));
    double phi_yy_01 = sec_der_dy(grid, func, grid.n_from_ij(i    , j+1));
    double phi_yy_11 = sec_der_dy(grid, func, grid.n_from_ij(i+1, j+1));
    double phi_yy = minmod(minmod(phi_yy_00, phi_yy_01),minmod(phi_yy_10, phi_yy_11));
    phi -= .5 * (y - y_j) * (y_jp1 - y) * phi_yy;

    return phi;
}


double bwd_dx(Grid2d & grid, std::vector<double> & func, int n){
    if ( int(func.size()) != (grid.get_N()*grid.get_M()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");
    else if ( grid.x_from_n(n) == grid.get_xmin() )
        return 0.;

    else
        return (func[n] - func[n-1])/grid.get_dx();

}
double fwd_dx(Grid2d & grid, std::vector<double> & func, int n){
    if ( int(func.size()) != (grid.get_N()*grid.get_M()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");
    else if ( grid.x_from_n(n) == grid.get_xmax() )
        return 0.;

    else
        return (func[n+1] - func[n])/grid.get_dx();

}
double bwd_dy(Grid2d & grid, std::vector<double> & func, int n){
    if ( int(func.size()) != (grid.get_N()*grid.get_M()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");
    else if ( grid.y_from_n(n) == grid.get_ymin() )
        return 0.;

    else
        return (func[n] - func[n-grid.get_N()])/grid.get_dy();

}
double fwd_dy(Grid2d & grid, std::vector<double> & func, int n){
    if ( int(func.size()) != (grid.get_N()*grid.get_M()) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");
    else if ( grid.y_from_n(n) == grid.get_ymax() )
        return 0.;

    else
        return (func[n+grid.get_N()] - func[n])/grid.get_dy();
}

double signum(double x){
    if ( x > 0. )
        return 1.;
    else if ( x < 0. )
        return -1.;
    else
        return 0.;
}

double ini_cond(double x, double y){
    return std::sqrt( pow((x-0.25), 2) + y*y ) - 0.2;
}

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