//
//  TrinomialTreePricer.cpp
//  TrinomialTreeII
//
//  Created by Changheng Chen on 9/10/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "TrinomialTreePricer.hpp"
#include <cmath>
#include <vector>

using namespace std;

// Return the max value of the two inputs
double max(double x, double y)
{
    return (x>y) ? x:y;
}

//====================================================================================================
// (1) Trinomial trees
double TT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    vector<vector<double>> S(2*N+1, vector<double>(N+1)); // Trinomial tree of S
    vector<vector<double>> V(2*N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                                   // Time interval
    double u = exp(v*sqrt(3*dt));                      // Factor of S moving up
    double d = exp(-v*sqrt(3*dt));                     // Factor of S moving down
    double pu = 1.0/6.0+(r-q-v*v/2)*sqrt(dt/(12*v*v)); // Risk neutral probabilities
    double pm = 2.0/3.0;
    double pd = 1.0/6.0-(r-q-v*v/2)*sqrt(dt/(12*v*v));
    
    // Create the Trinomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=2*j; i++)
        {
            if (i<j) S[i][j] = S0*pow(u,j-i);
            else if (i==j) S[i][j] = S0;
            else S[i][j] = S0*pow(d,i-j);
        }
    }
    
    // Calculate terminal payoffs
    for (i=0; i<=2*N; i++)
    {
        switch(PutCall)
        {
            case 'C': V[i][N] = max(S[i][N]-K, 0.0); break;
            case 'P': V[i][N] = max(K-S[i][N], 0.0); break;
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-1; j>=0; j--)
    {
        for (i=0; i<=2*j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = max(S[i][j]-K, exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    else
                        V[i][j] = max(K-S[i][j], exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    break;
                }
            }
        }
    }
    
    return V[0][0];
}

double TTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    vector<vector<double>> S(2*N+1, vector<double>(N+1)); // Trinomial tree
    vector<vector<double>> V(2*N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                                   // Time interval
    double u = exp(v*sqrt(3*dt));                      // Factor of S moving up
    double d = exp(-v*sqrt(3*dt));                     // Factor of S moving down
    double pu = 1.0/6.0+(r-q-v*v/2)*sqrt(dt/(12*v*v)); // Risk neutral probabilities
    double pm = 2.0/3.0;
    double pd = 1.0/6.0-(r-q-v*v/2)*sqrt(dt/(12*v*v));
    
    // Create the Trinomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=2*j; i++)
        {
            if (i<j) S[i][j] = S0*pow(u,j-i);
            else if (i==j) S[i][j] = S0;
            else S[i][j] = S0*pow(d,i-j);
        }
    }
    
    // Calculate terminal payoffs
    for (i=0; i<=2*N; i++)
    {
        switch(PutCall)
        {
            case 'C': V[i][N] = max(S[i][N]-K, 0.0); break;
            case 'P': V[i][N] = max(K-S[i][N], 0.0); break;
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-1; j>=0; j--)
    {
        for (i=0; i<=2*j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = max(S[i][j]-K, exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    else
                        V[i][j] = max(K-S[i][j], exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    break;
                }
            }
        }
    }

    // Compute the Greeks: [delta, gamma, theta]
    double G = numeric_limits<double>::quiet_NaN(); // Initialize with numeric_limits<double>::quiet_NaN()
    switch(Greek)
    {
        case 'D': G = (V[0][1]-V[2][1])/(S[0][1]-S[2][1]); break;
        //case 'G': G = (((V[0][2]-V[2][2])/(S[0][2]-S[2][2]))-((V[2][2]-V[4][2])/(S[2][2]-S[4][2])))/((S[0][2]-S[4][2])/2.0); break;
        case 'G': G = (((V[0][2]-V[2][2])/(S[0][2]-S[2][2]))-((V[2][2]-V[4][2])/(S[2][2]-S[4][2])))/(S[0][1]-S[2][1]); break;
        case 'T': G = (V[1][1]-V[0][0])/dt; break;
    }
    
    return G;
}

//====================================================================================================
// (2) Average Trinomial trees
double ATT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    return (TT(N, S0, K, T, q, r, v, PutCall, EuroAmer)+TT(N+1, S0, K, T, q, r, v, PutCall, EuroAmer))/2.0;
}

double ATTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    return (TTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek) + TTGreeks(N+1, S0, K, T, q, r, v, PutCall, EuroAmer, Greek))/2.0;
}

//====================================================================================================
// (3) TBS
double TBS(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    vector<vector<double>> S(2*N+1, vector<double>(N+1)); // Trinomial tree of S
    vector<vector<double>> V(2*N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                                   // Time interval
    double u = exp(v*sqrt(3*dt));                      // Factor of S moving up
    double d = exp(-v*sqrt(3*dt));                     // Factor of S moving down
    double pu = 1.0/6.0+(r-q-v*v/2)*sqrt(dt/(12*v*v)); // Risk neutral probabilities
    double pm = 2.0/3.0;
    double pd = 1.0/6.0-(r-q-v*v/2)*sqrt(dt/(12*v*v));
    
    // Create the Trinomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=2*j; i++)
        {
            if (i<j) S[i][j] = S0*pow(u,j-i);
            else if (i==j) S[i][j] = S0;
            else S[i][j] = S0*pow(d,i-j);
        }
    }
    
    // Calculate payoffs at (T-dt) using Black-Scholes formula
    for (i=0; i<=2*N-2; i++)
    {
        switch(EuroAmer)
        {
            case 'E': V[i][N-1] = BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall); break;
            case 'A':
            {
                if (PutCall=='C')
                    V[i][N-1] = max(S[i][N-1]-K, BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                else
                    V[i][N-1] = max(K-S[i][N-1], BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                break;
            }
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-2; j>=0; j--)
    {
        for (i=0; i<=2*j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = max(S[i][j]-K, exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    else
                        V[i][j] = max(K-S[i][j], exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    break;
                }
            }
        }
    }
    
    return V[0][0];
}

double TBSGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    vector<vector<double>> S(2*N+1, vector<double>(N+1)); // Trinomial tree of S
    vector<vector<double>> V(2*N+1, vector<double>(N+1)); // Option
    
    double dt = T/N;                                   // Time interval
    double u = exp(v*sqrt(3*dt));                      // Factor of S moving up
    double d = exp(-v*sqrt(3*dt));                     // Factor of S moving down
    double pu = 1.0/6.0+(r-q-v*v/2)*sqrt(dt/(12*v*v)); // Risk neutral probabilities
    double pm = 2.0/3.0;
    double pd = 1.0/6.0-(r-q-v*v/2)*sqrt(dt/(12*v*v));
    
    // Create the Trinomial price tree
    int i, j;                           // i, j--node counters on the vertical, horizontal axises
    for (j=0; j<=N; j++)
    {
        for(i=0; i<=2*j; i++)
        {
            if (i<j) S[i][j] = S0*pow(u,j-i);
            else if (i==j) S[i][j] = S0;
            else S[i][j] = S0*pow(d,i-j);
        }
    }
    
    // Calculate payoffs at (T-dt) using Black-Scholes formula
    for (i=0; i<=2*N-2; i++)
    {
        switch(EuroAmer)
        {
            case 'E': V[i][N-1] = BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall); break;
            case 'A':
            {
                if (PutCall=='C')
                    V[i][N-1] = max(S[i][N-1]-K, BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                else
                    V[i][N-1] = max(K-S[i][N-1], BlackScholes(S[i][N-1], K, dt, v, q, r, PutCall));
                break;
            }
        }
    }
    
    // Backward recursion through the tree to find option value
    for (j=N-2; j>=0; j--)
    {
        for (i=0; i<=2*j; i++)
        {
            switch(EuroAmer)
            {
                case 'E': V[i][j] = exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]); break;
                case 'A':
                {
                    if (PutCall=='C')
                        V[i][j] = max(S[i][j]-K, exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    else
                        V[i][j] = max(K-S[i][j], exp(-r*dt)*(pu*V[i][j+1]+pm*V[i+1][j+1]+pd*V[i+2][j+1]));
                    break;
                }
            }
        }
    }
    
    // Compute the Greeks: [delta, gamma, theta]
    double G = numeric_limits<double>::quiet_NaN(); // Initialize with numeric_limits<double>::quiet_NaN()
    switch(Greek)
    {
        case 'D': G = (V[0][1]-V[2][1])/(S[0][1]-S[2][1]); break;
        //case 'G': G = (((V[0][2]-V[2][2])/(S[0][2]-S[2][2]))-((V[2][2]-V[4][2])/(S[2][2]-S[4][2])))/((S[0][2]-S[4][2])/2.0); break;
        case 'G': G = (((V[0][2]-V[2][2])/(S[0][2]-S[2][2]))-((V[2][2]-V[4][2])/(S[2][2]-S[4][2])))/(S[0][1]-S[2][1]); break;    
        case 'T': G = (V[1][1]-V[0][0])/dt; break;
    }
    
    return G;
}

//====================================================================================================
// (4) TBSR
double TBSR(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer)
{
    
    if (N%2==1) N+=1; // Make sure N/2 is an integer
    return 2*TBS(N, S0, K, T, q, r, v, PutCall, EuroAmer)-TBS(N/2, S0, K, T, q, r, v, PutCall, EuroAmer);
}

double TBSRGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek)
{
    if (N%2==1) N+=1; // Make sure N/2 is an integer
    return 2*TBSGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek)-TBSGreeks(N/2, S0, K, T, q, r, v, PutCall, EuroAmer, Greek);
}

//====================================================================================================
// (0) Method selection function
double Pricer(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, int method)
{
    double price = numeric_limits<double>::quiet_NaN();
    switch(method)
    {
        case(1): price =   TT(N, S0, K, T, q, r, v, PutCall, EuroAmer); break;
        case(2): price =  ATT(N, S0, K, T, q, r, v, PutCall, EuroAmer); break;
        case(3): price =  TBS(N, S0, K, T, q, r, v, PutCall, EuroAmer); break;
        case(4): price = TBSR(N, S0, K, T, q, r, v, PutCall, EuroAmer); break;
    }
    return price;
}

double Greeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek, int method)
{
    double greek = numeric_limits<double>::quiet_NaN();
    switch(method)
    {
        case(1): greek =   TTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek); break;
        case(2): greek =  ATTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek); break;
        case(3): greek =  TBSGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek); break;
        case(4): greek = TBSRGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer, Greek); break;
    }
    return greek;
}
