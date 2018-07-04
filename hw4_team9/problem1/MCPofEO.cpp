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
#include "Statistics.hpp"
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

// (2.1) Monte Carlo pricing of european option
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

// (2.2.1) Pricing with Control Variate
double EO::MC_EO_CV(ULI N, char PutCall)
{
    std::vector<LDBL> z(N);         // Vector of random numbers
    std::vector<double> S(N), V(N); // Stock and option prices
    double b, VCV = .0;             // Coefficient and final option price
    
    // Generate N random numbers
    RNG MyRNG; z = MyRNG.BMM(N);
    
    // Calculate S and V
    for (ULI i=0; i<N; i++)
    {
        S[i] = data.S0*exp((data.r-data.q-data.v*data.v*0.5)*data.T+data.v*sqrt(data.T)*z[i]);
        switch(PutCall)
        {
            case 'C': V[i] = exp(-data.r*data.T)*max(S[i]-data.K, 0.0); break;
            case 'P': V[i] = exp(-data.r*data.T)*max(data.K-S[i], 0.0); break;
        }
    }
    
    // Calculate b
    b = covariance(S,V)/covariance(S,S);
    
    // Calculate final V
    for (ULI i=0; i<N; i++)
    {
        VCV += V[i]-b*(S[i]-exp(data.r*data.T)*data.S0);
    }
    
    return VCV/N;
}

// (2.2.2) Pricing with Antithetic Variables
double EO::MC_EO_AV(ULI N, char PutCall)
{
    std::vector<LDBL> z(N); // Vector of random numbers
    double S, VAV = .0;     // Stock and option price
    
    // Generate N random numbers using random numbers in [ 0, 1]
    RNG MyRNG; z = MyRNG.BMM(N);
    
    // Calculate S and V using z1 and z2
    for (ULI i=0; i<N; i++)
    {
        S = data.S0*exp((data.r-data.q-data.v*data.v*0.5)*data.T+data.v*sqrt(data.T)*z[i]);
        switch(PutCall)
        {
            case 'C': VAV += exp(-data.r*data.T)*max(S-data.K, 0.0); break;
            case 'P': VAV += exp(-data.r*data.T)*max(data.K-S, 0.0); break;
        }
        
        S = data.S0*exp((data.r-data.q-data.v*data.v*0.5)*data.T+data.v*sqrt(data.T)*(-z[i]));
        switch(PutCall)
        {
            case 'C': VAV += exp(-data.r*data.T)*max(S-data.K, 0.0); break;
            case 'P': VAV += exp(-data.r*data.T)*max(data.K-S, 0.0); break;
        }
    }
    
    return VAV/(2*N);
}

// (2.2.3) Pricing with Moment Matching
double EO::MC_EO_MM(ULI N, char PutCall)
{
    std::vector<LDBL> z(N);    // Vector of random numbers
    std::vector<double> S(N);  // Stock price
    double Si, Sbar, VMM = .0; // Intermediate, averaged stock prices and option price
    
    // Generate N random numbers
    RNG MyRNG; z = MyRNG.BMM(N);
    
    // Calculate S and its mean
    for (ULI i=0; i<N; i++)
    {
        S[i] = data.S0*exp((data.r-data.q-data.v*data.v*0.5)*data.T+data.v*sqrt(data.T)*z[i]);
    }
    Sbar = mean(S);
    
    // Update S and calculate V
    for (ULI i=0; i<N; i++)
    {
        Si = S[i]*(exp(data.r*data.T)*data.S0)/Sbar;
        switch(PutCall)
        {
            case 'C': VMM += exp(-data.r*data.T)*max(Si-data.K, 0.0); break;
            case 'P': VMM += exp(-data.r*data.T)*max(data.K-Si, 0.0); break;
        }
    }

    return VMM/N;
}

// (2.2.4) Pricing with Moment Matching and Control Variates
double EO::MC_EO_MMandCV(ULI N, char PutCall)
{
    std::vector<LDBL> z(N);         // Vector of random numbers
    std::vector<double> S(N), V(N); // Stock price, option price
    double Sbar;                    // Averaged stock prices
    double b = .0;                  // Coefficient b
    double VCVMM = .0;              // Option price

    // Generate N random numbers
    RNG MyRNG; z = MyRNG.BMM(N);
    
    // Calculate S
    for (ULI i=0; i<N; i++)
    {
        S[i] = data.S0*exp((data.r-data.q-data.v*data.v*0.5)*data.T+data.v*sqrt(data.T)*z[i]);
    }
    Sbar = mean(S);
    
    // Update S with averaged S, and calculate V
    for (ULI i=0; i<N; i++)
    {
        S[i] = S[i]*(exp(data.r*data.T)*data.S0)/Sbar;
        switch(PutCall)
        {
            case 'C': V[i] = exp(-data.r*data.T)*max(S[i]-data.K, 0.0); break;
            case 'P': V[i] = exp(-data.r*data.T)*max(data.K-S[i], 0.0); break;
        }
    }
    
    // Calculate b
    b = covariance(S,V)/covariance(S,S);
    
    // Calculate Vcv,mm
    for (ULI i=0; i<N; i++)
    {
        VCVMM += V[i]-b*(S[i]-exp(data.r*data.T)*data.S0);
    }
    
    return VCVMM/N;
}
