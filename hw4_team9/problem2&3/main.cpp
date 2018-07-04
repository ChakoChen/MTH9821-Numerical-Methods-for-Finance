//
//  main.cpp
//  PathDependentOption
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include "MCPofBO.hpp"
using namespace std;

double S01 = 25;
double S02 = 30;
double v1 = 0.3;
double v2 = 0.2;
double rho = 0.25;
double K = 50;
double T = 0.75;
double q = .0;
double r = 0.05;
char PutCall = 'C';

int main()
{
    //cout << "Please enter: 'C' for European call, 'P' for European put: ";
    //cin >> PutCall;
    cout << setprecision(12);

    double V; // Option price from Black-Scholes and Monte Carlo
    ULI N;    // Number of time steps of Monte Carlo simulations
    
    // (0) Creat option object
    BO MyOption; MyOption.SetBasket(S01, S02, v1, v2, rho, K, T, q, r);
    
    // (1) Monte Carlo pricing of basket option
    cout << "Basket Option: N, V(N) \n";
    for (ULI k=0; k<9; k++)
    {
        N = 10000*pow(2,k);
        V = MyOption.MC_BO(N, PutCall);
        cout << N << ", " << V << "\n";
    }
    cout << endl;
    
    // (2) Monte Carlo pricing of path-dependent basket option
    cout << "Path-Dependent Basket Option: n, m, V(n,m) \n";
    ULI n, m = 150;
    for (ULI k=0; k<10; k++)
    {
        n = (ULI)50*pow(2,k);
        V = MyOption.MC_Path_BO(n, m, PutCall);
        cout << n << ", " << m << ", " << V << "\n";
    }
    
    return 0;
}
