//
//  main.cpp
//  Heston Model
//
//  Created by Changheng Chen on 9/26/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include "MCPofHM.hpp"
using namespace std;

// Option parameters
double S0 = 50;     // Spot price
double K  = 50;     // Strike price
double T  = 0.50;   // Maturity
double v  = 0.30;   // Volatility
double q  = 0.00;   // Dividend
double r  = 0.05;   // Risk-free interest rate
char PutCall = 'P'; // Flag for put ('P') and call ('C') option

// Heston model parameters
double lambda = 3.0;     // Speed of mean-reversion
double vbar = 0.35*0.35; // Long term variance mean
double eta = 0.25;       // Standard deviation of the asset variance
double rho = -0.15;      // Correlation between the Wiener processes

int main()
{
    //cout << "Please enter: 'C' for European call, 'P' for European put: ";
    //cin >> PutCall;
    cout << setprecision(12);

    double V, sig;       // Option price and implied volatility
    ULI n, m = 175; // n: number of paths; m: number of time steps on each path
    
    // (0) Creat option object
    HM MyOption; MyOption.SetOption(S0, K, T, v, q, r);

    // (2) Monte Carlo pricing of path-dependent basket option
    cout << "Heston Model: n, m, V(m,n), sigma_BS \n";
    
    for (ULI k=0; k<6; k++)
    {
        n = (ULI)500*pow(2,k);
        V = MyOption.MC_HM(n, m, PutCall, lambda, vbar, eta, rho);
        sig = MyOption.implied_vol(V, PutCall);
        cout << n << ", " << m << ", " << V << ", " << sig << "\n";
    }
    
    return 0;
}
