//
//  MCPofEO.cpp
//  PathDependentOption: Monte Carlo Pricing of European Option
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <vector>
#include "MCPofEO.hpp"
#include "RandomNumberGenerator.hpp"

// Return the max value of the two inputs
double max(double x, double y)
{
    return (x>y) ? x:y;
}

// (0) Specification of option data
void EO::init()
{ // Initialize with arbitrary values
    
    data.S0 = 41;
    data.K  = 42;
    data.T  = 0.75;
    data.v  = 0.20;
    data.q  = 0.01;
    data.r  = 0.03;
}

void EO::SetOption(double S0, double K, double T, double v, double q, double r)
{
    data.S0 = S0;
    data.K  = K;
    data.T  = T;
    data.v  = v;
    data.q  = q;
    data.r  = r;
}

// (1) Constructors, destructor, and assignment operator
EO::EO(){init();};
EO::~EO(){};

// (2) Monte Carlo pricing of european option
std::tuple<double, double, double> EO::MC_EO(ULI N, char PutCall)
{ // N: number of paths; PutCall: flag for 'call' or 'put'
    
    std::vector<LDBL> z(N);           // Vector of random numbers
    RNG MyRNG; z = MyRNG.ITM(N);      // Generate N=n*m random numbers
    double S=0, V=0, delta=0, vega=0; // Underlying asset, option price, delta, and vega

    for (ULI i=0; i<N; i++)
    {
        S = data.S0*exp((data.r-data.q-data.v*data.v*0.5)*data.T+data.v*sqrt(data.T)*z[i]);
        switch(PutCall)
        {
            case 'C':
            {
                V += exp(-data.r*data.T)*max(S-data.K, 0.0);
                if (S-data.K>0)
                {
                    delta += exp(-data.r*data.T)*S/data.S0;
                    vega += S*exp(-data.r*data.T)*(-data.v*data.T+sqrt(data.T)*z[i]);
                }
                break;
            }
            case 'P':
            {
                V += exp(-data.r*data.T)*max(data.K-S, 0.0);
                if (S-data.K<0)
                {
                    delta += -exp(-data.r*data.T)*S/data.S0;
                    vega += -S*exp(-data.r*data.T)*(-data.v*data.T+sqrt(data.T)*z[i]);
                }
                break;
            }
        }
    }
    
    return std::make_tuple(V/N,delta/N,vega/N);
}
