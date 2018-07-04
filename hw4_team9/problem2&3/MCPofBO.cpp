//
//  MCPofBO.hpp
//  BasketOption: Monte Carlo Pricing of Basket Option
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <vector>
#include "MCPofBO.hpp"
#include "RandomNumberGenerator.hpp"

// Return the max value of the two inputs
double max(double x, double y)
{
    return (x>y) ? x:y;
}

// (0) Specification of option data
void BO::initBasket()
{ // Initialize with arbitrary values
    
    data.S01 = 25;
    data.S02 = 30;
    data.v1 = 0.3;
    data.v2 = 0.2;
    data.K = 50;
    data.T = 0.75;
    data.q = 0.0;
    data.r = 0.05;
    data.rho = 0.25;
}

void BO::SetBasket(double S01, double S02, double v1, double v2, double rho, double K, double T, double q, double r)
{
    data.S01 = S01;
    data.S02 = S02;
    data.v1 = v1;
    data.v2 = v2;
    data.rho = rho;
    data.K = K;
    data.T = T;
    data.q = q;
    data.r = r;
}

// (1) Constructors, destructor, and assignment operator
BO::BO(){initBasket();};
BO::~BO(){};

// (2.1) Monte Carlo pricing of basket option
double BO::MC_BO(ULI N, char PutCall)
{
    std::vector<LDBL> z(2*N); // Vector of random numbers
    double S = .0, V = .0;    // Underlying asset, option price, delta, and vega

    // Generate 2N random numbers
    RNG MyRNG; z = MyRNG.BMM(2*N);
    
    // Calculate S and V
    for (ULI i=0; i<N; i++)
    {
        S += data.S01*exp((data.r-data.q-data.v1*data.v1*0.5)*data.T + data.v1*sqrt(data.T)*z[2*i]);
        S += data.S02*exp((data.r-data.q-data.v2*data.v2*0.5)*data.T + data.v2*sqrt(data.T)*(data.rho*z[2*i]+sqrt(1-data.rho*data.rho)*z[2*i+1]));
        switch(PutCall)
        {
            case 'C': V += exp(-data.r*data.T)*max(S-data.K, 0.0); break;
            case 'P': V += exp(-data.r*data.T)*max(data.K-S, 0.0); break;
        }
        S = .0;
    }
    
    return V/N;
}

// (2.2) Monte Carlo pricing of path-dependent basket option
double BO::MC_Path_BO(ULI n, ULI m, char PutCall)
{ // n: number of paths; m: number of time steps on each path
    
    std::vector<LDBL> z(2*n*m); // Vector of random numbers
    double S1, S2, Smax;        // Record of stock prices and the max price along each path
    double V = .0;              // Option price
    
    // Generate N=n*m random numbers
    RNG MyRNG; z = MyRNG.BMM(2*n*m);
    
    // Calculate S1 and S2 along each path, and get the max to calculate V
    for (ULI i=0; i<n; i++)
    {
        S1 = data.S01; S2 = data.S02; Smax = S1+S2;
        for (ULI j=0; j<m; j++)
        {
            S1 = S1*exp((data.r-data.q-data.v1*data.v1*0.5)*(data.T/m) + data.v1*sqrt(data.T/m)*z[2*(i*m+j)]);
            S2 = S2*exp((data.r-data.q-data.v2*data.v2*0.5)*(data.T/m) + data.v2*sqrt(data.T/m)*(data.rho*z[2*(i*m+j)]+sqrt(1-data.rho*data.rho)*z[2*(i*m+j)+1]));
            Smax = max(Smax,S1+S2);
        }
        
        switch(PutCall)
        {
            case 'C': V += exp(-data.r*data.T)*max(Smax-data.K, 0.0); break;
            case 'P': V += exp(-data.r*data.T)*max(data.K-Smax, 0.0); break;
        }
    }

    return V/n;
}
