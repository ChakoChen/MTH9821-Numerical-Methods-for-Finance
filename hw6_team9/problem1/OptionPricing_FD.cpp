//
//  OptionPricing_FD.cpp
//  EuropeanOption_FD
//
//  Created by Changheng Chen on 10/15/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <algorithm>
#include <vector>
#include "BlackScholes.hpp"
#include "OptionPricing_FD.hpp"

using namespace std;


//========================================================================================================
// (0) Specification of option data
void FDM_OP::init()
{ // Initialize with arbitrary values
    
    data.S0 = 50;
    data.K  = 50;
    data.T  = 0.50;
    data.v  = 0.30;
    data.q  = 0.00;
    data.r  = 0.05;
    data.AmerEuro = 'E';
    data.PutCall = 'P';
}

void FDM_OP::SetOption(double S0, double K, double T, double v, double q, double r, char AmerEuro, char PutCall)
{
    data.S0 = S0;
    data.K  = K;
    data.T  = T;
    data.v  = v;
    data.q  = q;
    data.r  = r;
    data.AmerEuro = AmerEuro;
    data.PutCall = PutCall;
}


//========================================================================================================
// (1) Constructors, destructor, and assignment operator
FDM_OP::FDM_OP(){init();};
FDM_OP::~FDM_OP(){};


//========================================================================================================
// (2) Finite difference methods for option pricing
tuple<vector<vector<double>>, vector<double>, vector<double>> FDM_OP::FDM(double alpha_temp, int M)
{
    //----------------------------------------------------------------------------------------------------
    // (2.1) Calculate the parameters for the domain configuration
    double xleft  = log(data.S0/data.K)+(data.r-data.q-data.v*data.v*0.5)*data.T-3*data.v*sqrt(data.T);
    double xright = log(data.S0/data.K)+(data.r-data.q-data.v*data.v*0.5)*data.T+3*data.v*sqrt(data.T);
    double tfinal = data.T*data.v*data.v*0.5;
    
    double dt = tfinal/(double)M;
    int N = floor((xright-xleft)/sqrt(dt/alpha_temp));
    double dx = (xright-xleft)/(double)N;
    double alpha = dt/(dx*dx);
    
    double a = (data.r-data.q)/(data.v*data.v)-0.5;
    double b = pow((data.r-data.q)/(data.v*data.v)+0.5,2)+2*data.q/(data.v*data.v);

    vector<vector<double>> u(N+1, vector<double>(M+1)); // Vector to store solution (N+1 col, M+1 row)
    vector<double> statistics(3);                       // [error_pointwise, error_pointwise2, error_RMS]
    vector<double> greeks(3);                           // [delta, gamma, theta]
    
    //----------------------------------------------------------------------------------------------------
    // (2.2) Set the boundary conditions
    switch(data.PutCall)
    {
        case 'C':
        {
            for (int n=0; n<N+1; n++)
                u[n][0] = data.K*exp(a*(xleft+n*dx))*max(exp(xleft+n*dx)-1,0.);
            for (int m=0; m<M+1; m++)
            {
                u[0][m] = 0.;
                u[N][m] = data.K*exp(a*xright+b*(m*dt))*(exp(xright-2*data.q*(m*dt)/(data.v*data.v))-exp(-2*data.r*(m*dt)/(data.v*data.v)));
            }
            break;
        }
        case 'P':
        {
            for (int n=0; n<=N; n++)
                u[n][0] = data.K*exp(a*(xleft+n*dx))*max(1-exp(xleft+n*dx),0.);
            for (int m=0; m<=M; m++)
            {
                u[0][m] = data.K*exp(a*xleft+b*(m*dt))*(exp(-2*data.r*(m*dt)/(data.v*data.v))-exp(xleft-2*data.q*(m*dt)/(data.v*data.v)));
                u[N][m] = 0.;
            }
            break;
        }
    }
    
    //----------------------------------------------------------------------------------------------------
    // (2.3) Compute the values at the nodes of the inner domain
    for (int m=0; m<=M-1; m++)
    {
        for (int n=1; n<=N-1; n++)
        {
            u[n][m+1] = alpha*u[n-1][m] + (1-2*alpha)*u[n][m] + alpha*u[n+1][m];
        }
    }

    //----------------------------------------------------------------------------------------------------
    // (2.4) Calcuate option price and analytics
    //***** (2.4.0) Calculate Black-Scholes option price
    double Vexact = BlackScholes(data.S0, data.K, data.T, data.v, data.q, data.r, data.PutCall);
    
    //***** (2.4.1) First method for computing option price & corresponding pointwise relative error
    double xcompute = log(data.S0/data.K);
    int tag = 0;
    for (int n=0; n<N; n++)
    {
        if (xcompute>=(xleft+n*dx) && xcompute<(xleft+(n+1)*dx)) {tag = n; break;}
    }
    
    double x1 = xleft+tag*dx,   x2 = xleft+(tag+1)*dx;
    double S1 = data.K*exp(x1), S2 = data.K*exp(x2);  // S1: S_{i}; S2: S_{i+1}
    double V1 = exp(-a*x1-b*tfinal)*u[tag][M];        // V1: V_{i}
    double V2 = exp(-a*x2-b*tfinal)*u[tag+1][M];      // V2: V_{i+1}
    double Vapprox = ((S2-data.S0)*V1 + (data.S0-S1)*V2)/(S2-S1);
    
    double error_pointwise = abs(Vapprox-Vexact);
    
    //***** (2.4.2) Second method for computing option price & corresponding pointwise relative error
    double ucompute = ((x2-xcompute)*u[tag][M]+(xcompute-x1)*u[tag+1][M])/(x2-x1);
    double Vapprox2 = exp(-a*xcompute-b*tfinal)*ucompute;
    
    double error_pointwise2 = abs(Vapprox2-Vexact);
    
    //***** (2.4.3) Root-mean-squared error
    int NRMS = 0;
    double xk, Sk, Vapprox_k, Vexact_k, error_RMS = 0.;
    
    for (int n=0; n<=N; n++)
    {
        xk = xleft+n*dx;
        Sk = data.K*exp(xk);
        Vapprox_k = exp(-a*xk-b*tfinal)*u[n][M];
        Vexact_k = BlackScholes(Sk, data.K, data.T, data.v, data.q, data.r, data.PutCall);
        
        if (Vexact_k/data.S0>0.00001)
        {
            error_RMS += pow(Vapprox_k-Vexact_k, 2)/pow(Vexact_k, 2);
            NRMS += 1;
        }
    }
    error_RMS = sqrt(error_RMS/(double)NRMS);
    
    //----------------------------------------------------------------------------------------------------
    // (2.5) Calcuate the Greeks
    // Get variables for calculating Theta
    double V5 = exp(-a*x1-b*(tfinal-dt))*u[tag][M-1];   // V5: V_{i,delta_t}
    double V6 = exp(-a*x2-b*(tfinal-dt))*u[tag+1][M-1]; // V6: V_{i+1,delta_t}
    double Vapprox3 = ((S2-data.S0)*V5+(data.S0-S1)*V6)/(S2-S1);
    
    double theta = (Vapprox3-Vapprox)/(2*dt/pow(data.v,2));
    
    // Get variables for calculating Delta and Gamma
    x1 = xleft+(tag-1)*dx; x2 = xleft+(tag+2)*dx;
    double S3 = data.K*exp(x1), S4 = data.K*exp(x2); // S3: S_{i-1}; S4: S_{i+2}
    double V3 = exp(-a*x1-b*tfinal)*u[tag-1][M];     // V3: V_{i-1}
    double V4 = exp(-a*x2-b*tfinal)*u[tag+2][M];     // V4: V_{i+2}
    
    double delta = (V2-V1)/(S2-S1);
    double gamma = ((V4-V2)/(S4-S2)-(V1-V3)/(S1-S3))/((S4+S2)/2.-(S1+S3)/2.);
    
    //----------------------------------------------------------------------------------------------------
    // (2.6) Gather data to store in a tuple
    statistics = {error_pointwise, error_pointwise2, error_RMS};
    greeks = {delta, gamma, theta};
    
    return make_tuple(u, statistics, greeks);
}
