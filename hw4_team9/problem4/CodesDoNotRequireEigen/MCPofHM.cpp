//
//  MCPofHM.hpp
//  Heston Model: Monte Carlo Pricing Heston Model
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <vector>
#include "MCPofHM.hpp"
#include "BlackScholes.hpp"
#include "RandomNumberGenerator.hpp"

// Return the max value of the two inputs
double max(double x, double y)
{
    return (x>y) ? x:y;
}

// (0) Specification of option data
void HM::init()
{ // Initialize with arbitrary values
    
    data.S0 = 50;
    data.K  = 50;
    data.T  = 0.50;
    data.v  = 0.30;
    data.q  = 0.00;
    data.r  = 0.05;
}

void HM::SetOption(double S0, double K, double T, double v, double q, double r)
{
    data.S0 = S0;
    data.K  = K;
    data.T  = T;
    data.v  = v;
    data.q  = q;
    data.r  = r;
}

// (1) Constructors, destructor, and assignment operator
HM::HM(){init();};
HM::~HM(){};

// (2.1) Monte Carlo pricing for the Heston model
double HM::MC_HM(ULI n, ULI m, char PutCall, double lambda, double vbar, double eta, double rho)
{
    // n: number of paths
    // m: number of time steps on each path
    // PutCall: optin flag, 'P' for put and 'C' for call
    // lambda: teh speed of mean-reversion
    // vbar: long term variance mean
    // eta: standard deviation of the asset variance
    // rho: correlation of the Wiener processes

    std::vector<LDBL> z(2*n*m); // Vector of random numbers
    double var, S, V = .0;      // Variance, stock and option prices
    
    // Generate N=2*n*m random numbers
    RNG MyRNG; z = MyRNG.BMM(2*n*m);
    
    // Calculate S and V
    for (ULI i=0; i<n; i++)
    {
        S = data.S0; var = data.v*data.v;
        for (ULI j=0; j<m; j++)
        {
            S = S*exp((data.r-var*0.5)*(data.T/m)+sqrt(var*data.T/m)*z[2*(i*m+j)]);
            var = var-lambda*(var-vbar)*(data.T/m)+eta*sqrt(var*data.T/m)*(rho*z[2*(i*m+j)]+sqrt(1-rho*rho)*z[2*(i*m+j)+1]);
            var = max(var,.0);
        }
        
        switch(PutCall)
        {
            case 'C': V += exp(-data.r*data.T)*max(S-data.K, 0.0); break;
            case 'P': V += exp(-data.r*data.T)*max(data.K-S, 0.0); break;
        }
    }

    return V/n;
}

// (2.2) Black-Scholes implied volatility
double HM::implied_vol(double V, char PutCall)
{ // Newton's Method
    
    double X0 = 0.25; // Initial guess: 25% volatility
    double Xnew = X0, Xold = X0-1, tolerance = 1e-9; // Tolerance
    
    while (std::fabs(Xnew-Xold)>tolerance)
    {
        Xold = Xnew;
        Xnew = Xnew-(BlackScholes(data.S0, data.K, data.T, Xnew, data.q, data.r, PutCall)-V)/Vega(data.T, data.K, Xnew, data.r, data.q, data.S0);
    }
    
    return Xnew;
}
