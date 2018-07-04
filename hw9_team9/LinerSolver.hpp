//
//  LinerSolver.hpp
//  AmericanOption_FD
//
//  Created by Changheng Chen on 11/9/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef LinerSolver_hpp
#define LinerSolver_hpp

#include <vector>
#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

//a struct to combine the result of TriDiagonal LU decomposition
struct LUResult
{
    vector<vector<double>> Lmatrix;
    vector<vector<double>> Umatrix;
};

enum StoppingCriterion {consecutive, residual};

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


tuple<mat, mat> lu_no_pivoting(mat A)
{ // [L, U] = lu_no_pivoting(A)
    // A: nonsingular matrix with size n
    // L: lower triangular matrix with entries 1 on main diagonal
    // U: upper triangular matrix
    // A = LU
    
    long i, j, k;
    long n = A.rows();
    mat L(n,n), U(n,n);
    
    L.setZero(); U.setZero();
    for (i=0; i<=n-2; i++)
    {
        for (k=i; k<=n-1; k++)
        {
            U(i,k) = A(i,k);
            L(k,i) = A(k,i)/U(i,i);
        }
        for (j=i+1; j<=n-1; j++)
        {
            for (k=i+1; k<=n-1; k++)
            {
                A(j,k) -= L(j,i)*U(i,k);
            }
        }
    }
    L(n-1,n-1) = 1; U(n-1,n-1) = A(n-1,n-1);
    
    return make_tuple(L,U);
}

vec linear_solve_LU_no_pivoting(mat A, vec b)
{
    // A: nonsingular square matrix of size n with LU decomposition
    // b: column vector of size n
    // x = solution to Ax = b
    
    long n = A.rows();
    tuple<mat, mat> LU;
    mat L(n,n), U(n,n);
    vec x(n), y(n);

    LU = lu_no_pivoting(A);
    L = get<0>(LU);
    U = get<1>(LU);
    
    y = forward_subst(L, b);
    x = backward_subst(U, y);
    
    return x;
}


//To do Tridiagonal LU decomposition on A;
LUResult Tridiagonal_LU(vector<vector<double>> A) {
    long sz = A.size();
    vector<vector<double>> L(sz, vector<double>(sz, 0));
    vector<vector<double>> U(sz, vector<double>(sz, 0));
    
    for (int i = 0; i < sz-1; i++) {
        L[i][i] = 1;
        L[i + 1][i] = A[i + 1][i] / A[i][i];
        U[i][i] = A[i][i];
        U[i][i + 1] = A[i][i + 1];
        A[i + 1][i + 1] = A[i + 1][i + 1] - L[i + 1][i] * U[i][i + 1];
    }
    L[sz - 1][sz - 1] = 1;
    U[sz - 1][sz - 1] = A[sz - 1][sz - 1];
    
    LUResult Result = { L, U };
    return Result;
}


tuple <vec, uint> SOR(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion, const double & omega)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum;
    uint iter = 0;
    
    x.setZero();
    x0 = x_0;
    
    do{
        iter++;
        
        for (long i=0; i<n; i++)
        {
            sum = 0;
            for (long k=0; k<=i-1; k++)
            {
                sum += A(i,k)*x(k);
            }
            for (long k=i+1; k<=n-1; k++)
            {
                sum += A(i,k)*x0(k);
            }
            x(i) = max(x_0(i), (1-omega)*x0(i) + omega*(b(i)-sum)/A(i,i));
        }
        
        if (criterion == consecutive)
        {
            if ((x-x0).norm() < tolerance)
                break;
        }
        else
        {
            if ((b-A*x).norm() < tolerance)
                break;
        }
        
        x0 = x;
    } while (true);
    
    return make_tuple(x, iter);
}


#endif /* LinerSolver_hpp */

