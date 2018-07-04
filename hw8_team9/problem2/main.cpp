//
//  main.cpp
//  Exercise 2
//
//  Created by Changheng Chen on 10/30/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <Eigen/Dense>
#include <iomanip>
#include <vector>
#include <tuple>
#include <iostream>

using namespace Eigen;
using namespace std;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1, -1, uint> permutation;

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


tuple<permutation, mat, mat> lu_row_pivoting(mat A)
{ // [P, L, U] = lu_row_pivoting(A)
  // A: nonsingular matrix of size n
  // P: permutation matrix, stored as vector of its diagonal entries
  // L: lower triangular matrix with entries 1 on main diagonal
  // U: upper triangular matrix
  // PA = LU
    
    long i, j, k, m, l;
    long n = A.rows();
    permutation P(n);
    mat L(n,n), U(n,n);
    
    vector<double> vv, ww;
    long i_max;
    uint cc;
    double val;
    
    for (uint i=1; i<=n; i++){ P.indices()(i-1)=i; } // Initialize P = 1:n
    L.setIdentity();                                 // Initialize L as identity matrix
    U.setZero();                                     // Initialize U as matrix of 0s
    
    for (i=0; i<=n-2; i++)
    {
        // Find i_max, index of largest entry in |A(i:n,i)|
        i_max = i; val = abs(A(i,i));
        for (m=i+1; m<=n-1; m++)
        {
            if (abs(A(m,i))>val)
            {
                i_max = m;
                val = abs(A(m,i));
            }
        }
        
        // Switch rows i and i_max of A
        for (m=i; m<=n-1; m++)
        {
            vv.push_back(A(i,m));
            A(i,m) = A(i_max,m);
        }
        l = 0;
        for (m=i; m<=n-1; m++)
        {
            A(i_max, m) = vv[l]; l++;
        }
        
        // Update matrix P
        cc = P.indices()(i); P.indices()(i) = P.indices()(i_max); P.indices()(i_max) = cc;
        
        // switch rows i and i_max of L
        if (i>0)
        {
            for (m=0; m<=i-1; m++)
            {
                ww.push_back(L(i,m));
                L(i,m) = L(i_max,m);
            }
            l = 0;
            for (m=0; m<=i-1; m++)
            {
                L(i_max, m) = ww[l]; l++;
            }
        }
        
        for (j=i; j<=n-1; j++)
        {
            L(j,i) = A(j,i)/A(i,i); // compute column i of L
            U(i,j) = A(i,j);        // compute row i of U
        }
        for (j=i+1; j<=n-1; j++)
        {
            for (k=i+1; k<=n-1; k++)
            {
                A(j,k) = A(j,k)-L(j,i)*U(i,k);
            }
        }
        vv.clear(); ww.clear();
    }
    L(n-1,n-1) = 1; U(n-1,n-1) = A(n-1,n-1);
    
    return make_tuple(P, L, U);
}


int main()
{
    cout << setprecision(12);
    
    // (1) Test LU decomposition without row pivoting
    mat A(6,6);
    tuple<mat, mat> LU;
    A << -0.96, 0.82,-0.88, 0.61, 0.92, -0.3,
          -0.76,-1.17,-1.25,-0.56,-1.82,-0.05,
           0.73, 0.15, 0.67,-0.75, 0.10, 0.13,
          -0.17, 1.72, 0.55,-1.97,-2.07, 0.01,
          -0.39, 1.14,-0.34,-1.17, 0.64, 1.06,
          -0.09,-0.10,-0.34, 0.93, 1.60, 0.93;
    
    cout << "\n-----------------------------------------";
    cout << "\nTest LU decomposition without row pivoting\n";
    cout << "Here is the matrix A:\n" << A << endl;
    
    LU = lu_no_pivoting(A);
    cout << "Here is the matrix L:\n" << get<0>(LU) << endl;
    cout << "Here is the matrix U:\n" << get<1>(LU) << endl;
    
    
    // (2) Test LU decomposition with row pivoting
    tuple<permutation, mat, mat> PLU;
    A << -0.61,-0.15,-0.11, 1.99, 1.46,-1.51,
          0.35,-0.13,-1.48, 1.19, 0.10,-1.40,
          0.68, 1.49,-1.26, 0.92, 0.17, 0.43,
         -0.58,-0.29,-0.33,-1.10,-0.86,-0.09,
          1.53, 1.27,-1.18,-0.06, 2.33, 0.44,
          0.24,    0,-0.72,-0.09, 0.41,-0.49;
    
    cout << "\n-----------------------------------------";
    cout << "\nTest LU decomposition with row pivoting\n";
    cout << "Here is the matrix A:\n" << A << endl;
    
    PLU = lu_row_pivoting(A);
    cout << "Here is the matrix P:\n" << get<0>(PLU).indices() << endl;
    cout << "Here is the matrix L:\n" << get<1>(PLU) << endl;
    cout << "Here is the matrix U:\n" << get<2>(PLU) << endl;
    
    return 0;
}
