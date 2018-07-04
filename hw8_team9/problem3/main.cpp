//
//  main.cpp
//  Exercise 3
//
//  Created by Changheng Chen on 10/30/17.
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

//******************* (2) Cholesky decomposition *******************
mat cholesky(mat A)
{ // U = cholesky(A)
  // A: symmetric positive definite matrix of size n
  // U: upper triangular matrix s.t. (U^T)U=A
    
    long i, j, k;
    long n = A.rows();
    mat U(n,n);
    
    U.setZero();
    for (i=0; i<=n-2; i++)
    {
        U(i,i) = sqrt(A(i,i));
        for (k=i+1; k<=n-1; k++)
        {
            U(i,k) = A(i,k)/U(i,i);
        }
        for (j=i+1; j<=n-1; j++)
        {
            for (k=j; k<=n-1; k++)
            {
                A(j,k) -= U(i,j)*U(i,k);
            }
        }
    }
    U(n-1,n-1) = sqrt(A(n-1,n-1));
    
    return U;
}

//******** (3) Linear solver using Cholesky decomposition ********
vec liner_solve_cholesky(const mat & A, const vec & b)
{ // A: symmetric positive definite matrix of size n
  // b: column vector of size n
  // Ax = b
    
    long n = A.rows();
    mat U(n,n);
    vec x(n), y(n);

    U = cholesky(A);
    y = forward_subst(U.transpose(), b);
    x = backward_subst(U, y);
    
    return x;
}

int main()
{
    cout << setprecision(9);
    mat A(5,5), U(5,5);
    vec x(5), y(5), b(5);
    
    A <<  2.75, 1.61, 1.51, 1.22, 1.33,
          1.61, 3.93,-1.75, 2.27, 2.20,
          1.51,-1.75, 3.32,-1.10,-0.59,
          1.22, 2.27,-1.10, 2.58, 0.99,
          1.33, 2.20,-0.59, 0.99, 2.57;
    b << -0.02, 0.05,-0.45,-0.19,-1.84;
    
    cout << "\n-----------------------------------------";
    cout << "\nTest Cholesky decomposition\n";
    
    U = cholesky(A);
    x = liner_solve_cholesky(A, b);
    
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the matrix U:\n" << U << endl;
    cout << "Here is the vector x:\n" << x << endl;
    
    return 0;
}
