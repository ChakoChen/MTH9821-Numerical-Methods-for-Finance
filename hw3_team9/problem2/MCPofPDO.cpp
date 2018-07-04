//
//  MCPofPDO.cpp
//  PathDependentOption: Monte Carlo Pricing of a Path-Dependent Option
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <vector>
#include "MCPofPDO.hpp"
#include "BlackScholes.hpp"
#include "RandomNumberGenerator.hpp"

// Return the max value of the two inputs
double max(double x, double y)
{
    return (x>y) ? x:y;
}

// (0) Specification of option data
void PDO::init()
{ // Initialize with arbitrary values
    
    data.S0 = 39;
    data.K  = 39;
    data.B  = 35;
    data.T  = 0.75;
    data.v  = 0.25;
    data.q  = 0.01;
    data.r  = 0.02;
}

void PDO::SetOption(double S0, double K, double B, double T, double v, double q, double r)
{
    data.S0 = S0;
    data.K  = K;
    data.B  = B;
    data.T  = T;
    data.v  = v;
    data.q  = q;
    data.r  = r;
}

// (1) Constructors, destructor, and assignment operator
PDO::PDO(){init();};
PDO::~PDO(){};

// (2) Pricing of option
double PDO::exact_Cdao()
{
    return BlackScholes(data.S0, data.K, data.T, data.v, data.q, data.r, 'C') - pow(data.B/data.S0,2*((data.r-data.q)/(data.v*data.v))-0.5)*BlackScholes(data.B*data.B/data.S0, data.K, data.T, data.v, data.q, data.r, 'C');
}

double PDO::MC_Cdao(ULI n, ULI m)
{ // n: number of paths; m: number of time steps on each path
    
    std::vector<LDBL> z(n*m);      // Vector of random numbers
    double V = 0.0;                // Option price
    double S;                      // Record of underlying asset
    
    RNG MyRNG; z = MyRNG.ITM(n*m); // Generate N=n*m random numbers
    
    for (ULI i=0; i<n; i++)
    {
        S = data.S0;
        for (ULI j=0; j<m; j++)
        {
            S = S*exp((data.r-data.q-data.v*data.v*0.5)*(data.T/m)+data.v*sqrt(data.T/m)*z[i*m+j]);
            if (S <= data.B) continue; // If S <= B, V = 0; otherwise calculate V the following way
        }
        V += max(S-data.K,0.0);
    }

    return V*exp(-data.r*data.T)/n;
}
