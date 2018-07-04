//
//  BlackScholes.cpp
//  AmericanOption_FD
//
//  Created by Changheng Chen on 9/10/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "BlackScholes.hpp"
#include <cmath>
#include <math.h>

// PDF
double n(double x)
{
    return 1.0/sqrt(2.0 * M_PI) * exp(-x*x*0.5);
}

//Calculate N(t)
double N(double t)
{
    return 0.5*erfc(-t/sqrt(2));
}

// Put or call option from Black-Scholes model
double BlackScholes(double S, double K, double T, double v, double q, double r, char PutCall)
{
    double d1 = (log(S/K) + (r - q + (v*v)*0.5) * T)/(v * sqrt(T));
    double d2 = d1 - v * sqrt(T);
    
    return (PutCall=='C') ? ((S*exp(-q*T)*N(d1))-(K*exp(-r*T)*N(d2))) : ((K*exp(-r*T)*N(-d2))-(S*exp(-q*T)*N(-d1)));
}

// Delta of call and put
double CallDelta(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r-q + (v*v)*0.5) * T) / (v * sqrt(T));
    
    return exp(-q*T) * N(d1);
}

double PutDelta(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r-q + (v*v)*0.5) * T) / (v * sqrt(T));
    
    return -exp(-q*T) * N(-d1);
}

// Gamma
double Gamma(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r-q + (v*v)*0.5) * T) / (v * sqrt(T));
    
    return (n(d1) * exp(-q*T)) / (S* v * sqrt(T));
}

// Vega
double Vega(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r-q + (v*v)*0.5) * T) / (v * sqrt(T));
    
    
    return n(d1) * exp(-q*T) * sqrt(T) * S;
}

// Theta of call and put
double CallTheta(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r-q + (v*v)*0.5) * T) / (v * sqrt(T));
    double d2 = d1 - v * sqrt(T);
    double theta = -n(d1) * S*v*exp(-q*T)/(2*sqrt(T)) + q*S*exp(-q*T)*N(d1) - r*K*exp(-r*T)*N(d2);
    
    return theta;
}

double PutTheta(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r-q + (v*v)*0.5) * T) / (v * sqrt(T));
    double d2 = d1 - v * sqrt(T);
    double theta = -n(d1) * S*v*exp(-q*T)/(2*sqrt(T)) - q*S*exp(-q*T)*N(-d1) + r*K*exp(-r*T)*N(-d2);
    
    return theta;
}

double BlackScholesGreeks(double S, double K, double T, double v, double q, double r, char PutCall, char Greek)
{
    double G = std::numeric_limits<double>::quiet_NaN();
    
    switch(PutCall)
    {
        case('P'):
        {
            switch(Greek)
            {
                case('D'): G = PutDelta(T, K, v, r, q, S); break;
                case('G'): G = Gamma(T, K, v, r, q, S); break;
                case('T'): G = PutTheta(T, K, v, r, q, S); break;
                case('V'): G = Vega(T, K, v, r, q, S); break;
            }
            break;
        }
        case('C'):
        {
            switch(Greek)
            {
                case('D'): G = CallDelta(T, K, v, r, q, S); break;
                case('G'): G = Gamma(T, K, v, r, q, S); break;
                case('T'): G = CallTheta(T, K, v, r, q, S); break;
                case('V'): G = Vega(T, K, v, r, q, S); break;
            }
            break;
        }
    }
    
    return G;
}
