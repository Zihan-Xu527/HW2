//
// Created by Zihan Xu on 10/20/22.
//

#include "math_tools.h"
#include <cmath>

// suppose we have grid in 2D [xmin,xmax] x [ymin,ymax]
// find cell in which (x,y) belongs
// find weighted avg of values (?)

double bilinear_interpolation(Grid2d & grid, std::vector<double> & func, const double x, const double y){

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
        i = N;
    else if (x < xmin)
        i = 0;
    else
        i = floor((x - xmin)/dx);

    if (y > ymax)
        j = M;
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

double ENO_interpolation(Grid2d & grid, std::vector<double> & func, const double x, const double y){
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
        i = N;
    else if (x < xmin)
        i = 0;
    else
        i = floor((x - xmin)/dx);

    if (y > ymax)
        j = M;
    else if (y < ymin)
        j = 0;
    else
        j = floor((y - ymin)/dy);


//    int n = grid.n_from_ij(i,j);
//
//
//    std::cout << "i: " << i << " j: " << j << std::endl;
//    std::cout << "n: " << n << std::endl;

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

    double phi_xx_00 = sec_der_dx(grid, func, grid.n_from_ij(i    ,   j  ));
    double phi_xx_10 = sec_der_dx(grid, func, grid.n_from_ij(i+1,   j  ));
    double phi_xx_01 = sec_der_dx(grid, func, grid.n_from_ij(i    , j+1));
    double phi_xx_11 = sec_der_dx(grid, func, grid.n_from_ij(i+1, j+1));

    double phi_xx_minmod = minmod(minmod(phi_xx_00, phi_xx_01),minmod(phi_xx_10, phi_xx_11));
    phi -= .5 * (x - x_i) * (x_ip1 - x) * phi_xx_minmod;

    double phi_yy_00 = sec_der_dy(grid, func, grid.n_from_ij(i    ,   j  ));
    double phi_yy_10 = sec_der_dy(grid, func, grid.n_from_ij(i+1,   j  ));
    double phi_yy_01 = sec_der_dy(grid, func, grid.n_from_ij(i    , j+1));
    double phi_yy_11 = sec_der_dy(grid, func, grid.n_from_ij(i+1, j+1));

    double phi_yy_minmod = minmod(minmod(phi_yy_00, phi_yy_01),minmod(phi_yy_10, phi_yy_11));
    phi -= .5 * (y - y_j) * (y_jp1 - y) * phi_yy_minmod;

    return phi;

}
double MAX(double a, double b){
    if (a > b)
        return a;
    else
        return b;
}

int MAX(int a, int b){
    if (a > b)
        return a;
    else
        return b;
}

double minmod(double a, double b){
    if ( a*b < 0.0 )
        return 0.;
    else if ( std::abs(a) < std::abs(b) )
        return a;
    else
        return b;
}

double central_diff(double lo, double mid, double hi, double dx){
    return (lo - 2. * mid + hi) / (dx * dx);
}

double sec_der_dx(Grid2d & grid, std::vector<double> &func, int n)
{
    int N = grid.get_N();
    int M = grid.get_M();
    double dx = grid.get_dx();

    if ( int(func.size()) != (N*M) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

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

double sec_der_dy(Grid2d & grid, std::vector<double> &func, int n)
{
    int N = grid.get_N();
    int M = grid.get_M();
    double dy = grid.get_dy();

    if ( int(func.size()) != (N*M) )
        throw std::invalid_argument("ERROR: Dimension doesn't match!");

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