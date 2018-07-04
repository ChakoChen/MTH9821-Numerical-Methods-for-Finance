//
//  BinomialTreePricer.cpp
//  BinomialTreeII
//
//  Created by Changheng Chen on 9/3/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "BinomialTreePricer.hpp"
#include "BlackScholes.hpp"
#include <cmath>

using namespace std;

// Return the Max value of the two inputs
double Max(double x, double y)
{
    return (x>y) ? x:y;
}

//====================================================================================================
// (1) Binomial trees
double BT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    vector<vector<double>> S(N+1, vector<double>(N+1)); // Binomial tree of S
    vector<vector<double>> V(N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                    // Time interval
    double u = exp(v*sqrt(dt));         // Factor of S moving up
    double d = exp(-v*sqrt(dt));        // Factor of S moving down
    double p = (exp((r-q)*dt)-d)/(u-d); // Risk neutral probability
    
    // Create the binomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=j; i++)
        {
            S[i][j] = S0*pow(u,j-i)*pow(d,i);
        }
    }
    
    // Calculate terminal payoffs
    for (i=0; i<=N; i++)
    {
        switch(PutCall)
        {
            case 'C': V[i][N] = Max(S[i][N]-K, 0.0); break;
            case 'P': V[i][N] = Max(K-S[i][N], 0.0); break;
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-1; j>=0; j--)
    {
        for (i=0; i<=j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = Max(S[i][j]-K, exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    else
                        V[i][j] = Max(K-S[i][j], exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    break;
                }
            }
        }
    }
    
    return V[0][0];
}

double BTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    vector<vector<double>> S(N+1, vector<double>(N+1)); // Binomial tree of S
    vector<vector<double>> V(N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                    // Time interval
    double u = exp(v*sqrt(dt));         // Factor of S moving up
    double d = exp(-v*sqrt(dt));        // Factor of S moving down
    double p = (exp((r-q)*dt)-d)/(u-d); // Risk neutral probability
    
    // Create the binomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=j; i++)
        {
            S[i][j] = S0*pow(u,j-i)*pow(d,i);
        }
    }
    
    // Calculate terminal payoffs
    for (i=0; i<=N; i++)
    {
        switch(PutCall)
        {
            case 'C': V[i][N] = Max(S[i][N]-K, 0.0); break;
            case 'P': V[i][N] = Max(K-S[i][N], 0.0); break;
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-1; j>=0; j--)
    {
        for (i=0; i<=j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = Max(S[i][j]-K, exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    else
                        V[i][j] = Max(K-S[i][j], exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    break;
                }
            }
        }
    }

    // Compute the Greeks: [delta, gamma, theta]
    double G = numeric_limits<double>::quiet_NaN(); // Initialize with numeric_limits<double>::quiet_NaN()
    switch(Greek)
    {
        case 'D': G = (V[0][1]-V[1][1])/(S[0][1]-S[1][1]); break;
        case 'G': G = (((V[0][2]-V[1][2])/(S[0][2]-S[1][2]))-((V[1][2]-V[2][2])/(S[1][2]-S[2][2])))/((S[0][2]-S[2][2])/2.0); break;
        case 'T': G = (V[1][2]-V[0][0])/(2*dt); break;
    }
    
    return G;
}

//====================================================================================================
// (2) Average binomial trees
double ABT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    return (BT(N, S0, K, T, q, r, v, PutCall, EuroAmer)+BT(N+1, S0, K, T, q, r, v, PutCall, EuroAmer))/2.0;
}

double ABTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    return (BTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek) + BTGreeks(N+1, S0, K, T, q, r, v, PutCall, EuroAmer, Greek))/2.0;
}

//====================================================================================================
// (3) BBS
double BBS(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    vector<vector<double>> S(N+1, vector<double>(N+1)); // Binomial tree of S
    vector<vector<double>> V(N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                    // Time interval
    double u = exp(v*sqrt(dt));         // Factor of S moving up
    double d = exp(-v*sqrt(dt));        // Factor of S moving down
    double p = (exp((r-q)*dt)-d)/(u-d); // Risk neutral probability
    
    // Create the binomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=j; i++)
        {
            S[i][j] = S0*pow(u,j-i)*pow(d,i);
        }
    }
    
    // Calculate payoffs at (T-dt) using Black-Scholes formula
    for (i=0; i<=N-1; i++)
    {
        switch(EuroAmer)
        {
            case 'E': V[i][N-1] = BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall); break;
            case 'A':
            {
                if (PutCall=='C')
                    V[i][N-1] = Max(S[i][N-1]-K, BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                else
                    V[i][N-1] = Max(K-S[i][N-1], BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                break;
            }
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-2; j>=0; j--)
    {
        for (i=0; i<=j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = Max(S[i][j]-K, exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    else
                        V[i][j] = Max(K-S[i][j], exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    break;
                }
            }
        }
    }
    
    return V[0][0];
}

double BBSGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    vector<vector<double>> S(N+1, vector<double>(N+1)); // Binomial tree of S
    vector<vector<double>> V(N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                    // Time interval
    double u = exp(v*sqrt(dt));         // Factor of S moving up
    double d = exp(-v*sqrt(dt));        // Factor of S moving down
    double p = (exp((r-q)*dt)-d)/(u-d); // Risk neutral probability
    
    // Create the binomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=j; i++)
        {
            S[i][j] = S0*pow(u,j-i)*pow(d,i);
        }
    }
    
    // Calculate payoffs at (T-dt) using Black-Scholes formula
    for (i=0; i<=N-1; i++)
    {
        switch(EuroAmer)
        {
            case 'E': V[i][N-1] = BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall); break;
            case 'A':
            {
                if (PutCall=='C')
                    V[i][N-1] = Max(S[i][N-1]-K, BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                else
                    V[i][N-1] = Max(K-S[i][N-1], BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                break;
            }
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-2; j>=0; j--)
    {
        for (i=0; i<=j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = Max(S[i][j]-K, exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    else
                        V[i][j] = Max(K-S[i][j], exp(-r*dt)*(p*V[i][j+1]+(1-p)*V[i+1][j+1]));
                    break;
                }
            }
        }
    }
    
    // Compute the Greeks: [delta, gamma, theta]
    double G = numeric_limits<double>::quiet_NaN(); // Initialize with numeric_limits<double>::quiet_NaN()
    switch(Greek)
    {
        case 'D': G = (V[0][1]-V[1][1])/(S[0][1]-S[1][1]); break;
        case 'G': G = (((V[0][2]-V[1][2])/(S[0][2]-S[1][2]))-((V[1][2]-V[2][2])/(S[1][2]-S[2][2])))/((S[0][2]-S[2][2])/2.0); break;
        case 'T': G = (V[1][2]-V[0][0])/(2*dt); break;
    }
    
    return G;
}

//====================================================================================================
// (4) BBSR
double BBSR(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    if (N%2==1) N+=1; // Make sure N/2 is an integer
    return 2*BBS(N, S0, K, T, q, r, v, PutCall, EuroAmer)-BBS(N/2, S0, K, T, q, r, v, PutCall, EuroAmer);
}

double BBSRGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    if (N%2==1) N+=1; // Make sure N/2 is an integer
    return 2*BBSGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek)-BBSGreeks(N/2, S0, K, T, q, r, v, PutCall, EuroAmer, Greek);
}
