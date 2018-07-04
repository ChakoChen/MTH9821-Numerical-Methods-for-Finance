//
//  main.cpp
//  Exercise 1
//
//  Created by Changheng Chen on 10/29/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1, -1, uint> permutation;

//******************* (1.1) Forward substitution *******************
vec forward_subst(const mat & L, const vec & b)
{ /* Solve x for Lx = b;
     L: nonsingular lower triangular matrix of size n
     b: column vector of size n
   */
    
    long n = L.rows();
    vec x(n);
    x.setZero();
    double sum;
    
    if (L.rows() != L.cols() || L.cols() != b.size())
    {
        cout << "forward_subst: dimension mismatch\n";
        return x;
    }

    x(0) = b(0) / L(0, 0);
    for (long j=1; j<=n-1; j++)
    {
        sum = 0.;
        for (long k=0; k<j; k++)
        {
            sum += L(j, k) * x(k);
        }
        x(j) = (b(j)-sum)/L(j, j);
    }
    
    return x;
}

//******************* (1.2) Backward substitution *******************
vec backward_subst(const mat & U, const vec & b)
{ /* Solve x for Ux = b;
   U: nonsingular upper triangular matrix of size n
   b: column vector of size n
   */
    
    long n = U.rows();
    vec x(n);
    x.setZero();
    double sum;
    
    if (U.rows() != U.cols() || U.cols() != b.size())
    {
        cout << "backward_subst: dimension mismatch\n";
        return x;
    }
    
    x(n-1) = b(n-1)/U(n-1, n-1);
    for (long j=n-2; j>=0; j--)
    {
        sum = 0.;
        for (long k=j+1; k<=n-1; k++)
        {
            sum += U(j,k)*x(k);
        }
        x(j) = (b(j)-sum)/U(j,j);
    }
    
    return x;
}

//******************* (2.2) Forward substitution (band) *******************
vec forward_subst_band(const mat & L, const uint bands, const vec & b)
{ /* Solve x for Lx = b;
   L: nonsingular lower triangular banded matrix of band 'bands' of size n
   b: column vector of size n
   */
    
    long n = L.rows();
    vec x(n);
    x.setZero();
    
    if (L.rows() != L.cols() || L.cols() != b.size() || bands>L.rows())
    {
        cout << "forward_subst: dimension mismatch\n";
        return x;
    }
    
    x = b;
    for (long j=0; j<=n-1; j++)
    {
        x(j) /= L(j,j);
        for (long k=j+1; k<std::min(j+1+bands, n); k++)
        {
            x(k) -= L(k,j)*x(j);
        }
    }
    
    return x;
}

//******************* (2.2) Backward substitution (band) *******************
vec backward_subst_band(const mat & U, const uint bands, const vec & b)
{ /* Solve x for Ux = b;
   U: nonsingular upper triangular banded matrix of band 'bands' of size n
   b: column vector of size n
   */
    
    long n = U.rows();
    vec x(n);
    x.setZero();
    
    if (U.rows() != U.cols() || U.cols() != b.size() || bands>U.rows())
    {
        cout << "forward_subst: dimension mismatch\n";
        return x;
    }
    
    x = b;
    for (long j=n-1; j>=0; j--)
    {
        x(j) /= U(j, j);
        for (long k=std::max((long)0, j-bands); k<j; k++)
        {
            x(k) -= U(k,j) * x(j);
        }
    }
    
    return x;
}



int main()
{
    cout << setprecision(9);
    
    // (1.1) Test forward_subst(const mat & L, const vec & b)
    mat A1(7,7);
    vec b1(7);
    
    A1 << -1.73,    0,    0,    0,    0,    0,   0,
          -0.01, 0.21,    0,    0,    0,    0,   0,
          -0.28, 0.76, 0.62,    0,    0,    0,   0,
           1.02, 0.74,-0.62,-0.14,    0,    0,   0,
          -0.57, 0.09,-1.40,-0.69,-0.68,    0,   0,
           0.38,-0.28,-0.81,-0.27,-0.44,-0.02,   0,
          -1.40,-1.38,-0.14,-2.60, 1.54, 1.65,0.42;
    b1 <<  3.63, 0.42, 2.15,   -1, 0.46, -1.7,-3.24;
    
    cout << "\n-----------------------------------------------\n"
         << "Test forward substitution \n";
    cout << "Here is the vector b:\n" << b1 << endl;
    cout << "Here is the matrix A:\n" << A1 << endl;
    vec x11 = A1.colPivHouseholderQr().solve(b1);
    vec x12 = forward_subst(A1, b1);
    cout << "The solution using Eigen intrinsic function:\n" << x11 << endl;
    cout << "The solution using forward_subst(L, b):\n" << x12 << endl;
    

    // (1.2) Test backward_subst(const mat & L, const vec & b)
    mat A2(5,5);
    vec b2(5);
    
    A2 <<  1.3, -4.9,  6.6,   2.4,-11.3,
             0, -2.7,-14.7, -13.5,  2.2,
             0,    0, 19.1,   0.8, -9.1,
             0,    0,    0,   0.5, 10.5,
             0,    0,    0,     0, 15.4;
    b2 << 1.56,-0.71, -1.1, -1.65, -0.36;
    
    cout << "\n-----------------------------------------------\n"
         << "Test backward substitution \n";
    cout << "Here is the vector b:\n" << b2 << endl;
    cout << "Here is the matrix A:\n" << A2 << endl;
    vec x21 = A2.colPivHouseholderQr().solve(b2);
    vec x22 = backward_subst(A2, b2);
    cout << "The solution using Eigen intrinsic function:\n" << x21 << endl;
    cout << "The solution using backward_subst(L, b):\n" << x22 << endl;
    
    
    // (2.2) Test forward_subst_band(const mat & L, const uint bands, const vec & b)
    mat A(3,3);
    vec b(3);
    
    A << 1,0,0, 4,5,0, 0,8,10;
    b << 3,3,4;
    
    cout << "\n-----------------------------------------------\n"
         << "Test forward substitution for banded matrix\n";
    cout << "Here is the vector b:\n" << b << endl;
    cout << "Here is the matrix A:\n" << A << endl;
    uint bands = 2;
    vec x1 = A.colPivHouseholderQr().solve(b);
    vec x2 = forward_subst_band(A, bands, b);
    cout << "The solution using Eigen intrinsic function:\n" << x1 << endl;
    cout << "The solution using forward_subst_band(L, band, b):\n" << x2 << endl;
    
    
    // (2.2) Test backward_subst_band(const mat & L, const uint bands, const vec & b)
    A << 1,2,0, 0,5,6, 0,0,10;
    
    cout << "\n-----------------------------------------------\n"
         << "Test backward substitution for banded matrix\n";
    cout << "Here is the vector b:\n" << b << endl;
    cout << "Here is the matrix A:\n" << A << endl;
    bands = 2;
    x1 = A.colPivHouseholderQr().solve(b);
    x2 = backward_subst_band(A, bands, b);
    cout << "The solution using Eigen intrinsic function:\n" << x1 << endl;
    cout << "The solution using backward_subst_band(L, band, b):\n" << x2 << endl;
    
    return 0;
}
