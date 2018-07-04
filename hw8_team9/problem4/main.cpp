//
//  main.cpp
//  Exercise 4
//
//  Created by Changheng Chen on 10/31/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <vector>
#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1, -1, uint> permutation;

enum StoppingCriterion {consecutive, residual};

// (1) Jacobi iterative methods: with and without initial guess x0 being given
tuple<vec, uint> jacobi(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum;
    uint k = 1;
    
    x = x_0; x0 = x_0;

    do {
        for (long i=0; i<n; i++)
        {
            sum = 0;
            for (long j=0; j<n; j++)
            {
                if (i != j)
                    sum += -A(i,j)*x0[j];
            }
            x[i] = (b[i]+sum)/A(i,i);
        }
        
        k++;
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
    
    return make_tuple(x, k);
}

tuple<vec, uint> jacobi(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum;
    uint k = 0;
    
    x.setZero();
    x0 = x;

    do {
        for (long i=0; i<n; i++)
        {
            sum = 0;
            for (long j=0; j<n; j++)
            {
                if (i != j)
                    sum += -A(i,j)*x0[j];
            }
            x[i] = (b[i]+sum)/A(i,i);
        }
        
        k++;
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
    
    cout << "The norm is " << (x-x0).norm() << "\n";
    return make_tuple(x, k);
}

// (2) Gauss-Siedel method with and without initial guess x0 being given
tuple <vec, uint> gauss_seidel(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum;
    uint iter = 0;
    
    x.setZero();
    x0 = x;
    
    do {
        iter++;
        
        for (long i = 0; i < n; i++)
        {
            sum = 0;
            for (long k = 0; k <= i-1; k++)
            {
                sum += A(i,k)*x(k);
            }
            for (long k=i+1; k<=n-1; k++)
            {
                sum += A(i,k)*x0(k);
            }
            x[i] = (b[i]-sum)/A(i,i);
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
    
    cout << "The norm is " << (x-x0).norm() << "\n";
    
    return make_tuple(x, iter);
}


tuple <vec, uint> gauss_seidel(const mat & A, const vec & b, const vec & x_0 , const double tolerance, const StoppingCriterion criterion)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum;
    uint iter = 0;
    
    x.setZero();
    x = x_0;
    x0 = x_0;
    
    do {
        iter++;
        
        for (long i = 0; i < n; i++)
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
            x[i] = (b[i]-sum)/A(i,i);
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

// (3) SOR method with and without initial guess being given
tuple <vec, uint> SOR(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion, const double & omega)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum = 0;
    uint iter = 0;
    
    x.setZero();
    x0 = x;
    
    do {
        iter++;
        
        for (long i = 0; i < n; i++)
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
            x[i] = (1-omega)*x0[i]+omega*(b[i]-sum)/A(i,i);
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
    
    cout << "The norm is " << (x-x0).norm() << "\n";
    
    return make_tuple(x, iter);
}

tuple <vec, uint> SOR(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion, const double & omega)
{
    long n = b.size();
    vec x(n), x0(n);
    double sum;
    uint iter = 0;
    
    x.setZero();
    x = x_0;
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
            x[i] = (b[i]-sum)/A(i,i);
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


int main()
{
    cout << setprecision(9);
    
    cout << "========================================\nTest Jacobi Method: \n";
    mat A(4, 4);
    vec b(4);
    
    A << 5,-2, 1, 1,
         2, 7, 2,-1,
         1, 2, 7, 1,
        -1, 1, 2, 8;
    b << 16, 11, 16, 19;
    
    tuple<vec, uint> res1 = jacobi(A, b, 0.000001, consecutive);
    cout << "Iteration times: " << get<1>(res1)
         << "\nSolution to the linear system: \n" << std::get<0>(res1);
    
    cout << "\n\n========================================\nTest Gauss-Siedel Method: \n";
    A.resize(5, 5);
    b.resize(5);
    
    A << 10,  3, -2,  1,  1,
          2,  9,  2,  2, -1,
          1,  2, 13,  6,  1,
         -1,  1,  2, 24,  5,
          3,  3,  2, -4, 19;
    b << 16, 11, 16, 19, 20;
    
    tuple<vec, uint> res2 = gauss_seidel(A, b, 0.000001, consecutive);
    cout << "Iteration times: " << get<1>(res2)
         << "\nSolution to the linear system: \n" << std::get<0>(res2);
    
    cout << "\n\n========================================\nTest SOR Method: \n";
    A.resize(6, 6);
    b.resize(6);

    A << 15,  5,  3, -2,  1,  1,
          2, 26,  6,  2,  2, -1,
          1,  2, 18,  3,  6,  1,
         -1,  1,  2, 28,  5,  5,
          3,  3,  2, -4, 36,  8,
         -9,  8,  3,  1,  2, 45;
    b << 16, 11, 16, 19, 20,-19;
    
    tuple<vec, uint> res3 = SOR(A, b, 0.000001, consecutive, 1.05);
    cout << "Iteration times: " << get<1>(res3)
         << "\nSolution to the linear system: \n" << std::get<0>(res3) << endl;

    return 0;
}
