//
//  OptionPricing_FD.cpp
//  AmericanOption_FD
//
//  Created by Changheng Chen on 10/15/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
#include "BlackScholes.hpp"
#include "LinerSolver.hpp"
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
    data.AmerEuro = 'A';
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
// (1) Constructors and destructor
FDM_OP::FDM_OP(){init();};
FDM_OP::~FDM_OP(){};


//========================================================================================================
// (2) Finite difference methods for option pricing
tuple<vector<vector<double>>, vector<double>, vector<double>> FDM_OP::FDM(double alpha_temp, double omega, int M)
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

    vector<vector<double>> early_ex_premium(N-1, vector<double>(M)); // Early exercise premium matrix
    vector<vector<double>> u(N+1, vector<double>(M+1)); // Matrix to store solution (N+1 col, M+1 row)
    vector<double> statistics(4);  // [error_pointwise, error_pointwise2, VarRed, error_pointwise_VarRed]
    vector<double> greeks(3);      // [delta, gamma, theta]
    
    //----------------------------------------------------------------------------------------------------
    // (2.2) Set the initial and boundary conditions
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
                if (data.AmerEuro=='E')
                    u[0][m] = data.K*exp(a*xleft+b*m*dt)*(exp(-2*data.r*(m*dt)/(data.v*data.v))-exp(xleft-2*data.q*m*dt/(data.v*data.v)));
                else
                    u[0][m] = data.K*exp(a*xleft+b*(m*dt))*(1-exp(xleft));
                
                u[N][m] = 0.;
            }
            break;
        }
    }
    
    //----------------------------------------------------------------------------------------------------
    // (2.3) Compute early exercise premium in the inner domain
    for (int m=0; m<=M-1; m++)
    {
        for (int n=0; n<=N-2; n++)
            early_ex_premium[n][m] = data.K*exp(a*(xleft+(n+1)*dx)+b*(m+1)*dt) * max(1-exp(xleft+(n+1)*dx),0.);
    }

    //----------------------------------------------------------------------------------------------------
    // (2.4) Compute values at the nodes of inner domain using Crank-Nicolson scheme with SOR solver
    mat A(N-1, N-1);
    vec b_vec(N-1);
    vec x_0(N-1);
    
    tuple <vec, uint> res;
    double tolerance = 1e-6;
    
    // initialize A
    A.setZero();
    for (int n=0; n<=N-2; n++)
    {
        if (n==0)
        {
            A(n,n)   = 1+alpha;
            A(n,n+1) = -alpha/2;
        }
        else if(n==N-2)
        {
            A(n,n-1) = -alpha/2;
            A(n,n)   = 1+alpha;
        }
        else
        {
            A(n,n-1) = -alpha/2;
            A(n,n)   = 1+alpha;
            A(n,n+1) = -alpha/2;
        }
    }
    
    for (int m=0; m<=M-1; m++)
    {
        // initialize b_vec and initial guess x_0
        for (int n=0; n<=N-2; n++)
        {
            b_vec(n) = alpha/2.0 * u[n+2][m] + (1-alpha)*u[n+1][m] + alpha/2.0*u[n][m];
            x_0(n) = early_ex_premium[n][m];
        }
        
        // compute u time-step by time-step
        res = SOR(A, b_vec, x_0, tolerance, consecutive, omega);
        
        x_0 = get<0>(res);
        for (int n=0; n<=N-2; n++)
            u[n+1][m+1] = x_0(n);
    }
    
    //----------------------------------------------------------------------------------------------------
    // (2.5) Calcuate option price and analytics
    //***** (2.4.0) Calculate the "exact" option price from binomial tree mthod
    double Vexact = 3.3045362802172642;
    
    //***** (2.5.1) First method for computing option price & corresponding pointwise relative error
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
    
    //***** (2.5.2) Second method for computing option price & corresponding pointwise relative error
    double ucompute = ((x2-xcompute)*u[tag][M]+(xcompute-x1)*u[tag+1][M])/(x2-x1);
    double Vapprox2 = exp(-a*xcompute-b*tfinal)*ucompute;
    
    double error_pointwise2 = abs(Vapprox2-Vexact);
    
    //----------------------------------------------------------------------------------------------------
    // (2.6) Calcuate the Greeks
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
    // (2.7) Calcuate variance reduction values
    double VBS = BlackScholes(data.S0, data.K, data.T, data.v, data.q, data.r, data.PutCall);
    double Vapprox_Eur = FDM(alpha_temp, M);
    double VarRed = Vapprox + (VBS-Vapprox_Eur);
    double error_pointwise_VarRed = abs(VarRed-Vexact);
    
    //----------------------------------------------------------------------------------------------------
    // (2.8) Gather data to store in a tuple
    statistics = {error_pointwise, error_pointwise2, VarRed, error_pointwise_VarRed};
    greeks = {delta, gamma, theta};
    
    return make_tuple(u, statistics, greeks);
}


//return be values calculated by CrankNicolson Method
double FDM_OP::FDM(double alpha_temp, int M)
{
    double xleft  = log(data.S0/data.K)+(data.r-data.q-data.v*data.v*0.5)*data.T-3*data.v*sqrt(data.T);
    double xright = log(data.S0/data.K)+(data.r-data.q-data.v*data.v*0.5)*data.T+3*data.v*sqrt(data.T);
    double tfinal = data.T*data.v*data.v*0.5;
    
    double dt = tfinal/(double)M;
    int N = floor((xright-xleft)/sqrt(dt/alpha_temp));
    double dx = (xright-xleft)/(double)N;
    double alpha = dt/(dx*dx);
    
    double a = (data.r-data.q)/(data.v*data.v)-0.5;
    double b = pow((data.r-data.q)/(data.v*data.v)+0.5,2)+2*data.q/(data.v*data.v);
    
    //construct maxtrix A(size: N-1)
    vector<vector<double>>  A(N - 1, vector<double>(N - 1, 0));
    for (int i = 1; i < N - 2; i++) {
        A[i][i] = 1 + alpha;
        A[i][i - 1] = -0.5 * alpha;
        A[i][i + 1] = -0.5 * alpha;
    }
    A[0][0] = 1 + alpha;
    A[0][1] = -0.5 * alpha;
    A[N - 2][N - 2] = 1 + alpha;
    A[N - 2][N - 3] = -0.5 * alpha;
    
    //Get Tridiagonal_LU decomposition of A
    LUResult LandU = Tridiagonal_LU(A);
    vector<vector<double>> L = LandU.Lmatrix;
    vector<vector<double>> U = LandU.Umatrix;
    
    //construct the result for each point
    vector<vector<double>> Result(M + 1, vector<double>(N + 1));
    
    //calculte for tao = 0 (m=0)
    for (int i = 0; i < N + 1; i++)
    {
        double x_i = xleft + i*dx;
        Result[0][i] = data.K*exp(a*x_i)*max(1 - exp(x_i), 0.0);
    }
    //calculate for x_left(n=0) and x_right(n=N)
    for (int i = 0; i < M + 1; i++)
    {
        double tao_i = i*dt;
        //x_left
        Result[i][0] = data.K*exp(a*xleft + b*tao_i)*(exp(-2 * tao_i*data.r / pow(data.v, 2)) - exp(xleft - 2 * data.q*tao_i / pow(data.v, 2)));
        //x_tight
        Result[i][N] = 0;
    }
    
    //calculate for m=1:M
    for (int m = 1; m < M + 1; m++)
    {
        //construct vector b(size:N-1);
        vector<double> b(N - 1);
        for (int i = 1; i < N - 2; i++) {
            b[i] = 0.5*alpha*Result[m - 1][i + 2] + (1 - alpha)*Result[m - 1][i + 1] + 0.5*alpha*Result[m - 1][i];
        }
        b[0] = 0.5*alpha*Result[m - 1][2] + (1 - alpha)*Result[m - 1][1]
        + 0.5*alpha*Result[m - 1][0] + 0.5*alpha*Result[m][0];
        b[N - 2] = 0.5*alpha*Result[m - 1][N] + (1 - alpha)*Result[m - 1][N - 1]
        + 0.5*alpha*Result[m - 1][N - 2] + 0.5*alpha*Result[m][N];
        
        //Solve for LUx=b : find y s.t. Ly=b (size of y: N-1)
        vector<double> Y(N - 1);
        Y[0] = b[0];
        for (int i = 1; i < N - 1; i++) {
            Y[i] = b[i] - L[i][i - 1] * Y[i - 1];
        }
        
        //Solve for Ux=y £º find x s.t. Ux= y(size of x: N-1)
        Result[m][N - 1] = Y[N - 2] / U[N - 2][N - 2];
        for (int j = N - 3; j >= 0; j--) {
            Result[m][j + 1] = (Y[j] - U[j][j + 1] * Result[m][j + 2]) / U[j][j];
        }
    }
    //Analysis Convergence
    double x_comp = log(data.S0 / data.K);
    int i = floor((x_comp - xleft) / dx);
    double x_1 = i*dx + xleft;
    double x_2 = (i + 1)*dx + xleft;
    
    double S_1 = data.K*exp(x_1);
    double S_2 = data.K*exp(x_2);
    double V_1 = exp(-a*x_1 - b*tfinal)*Result[M][i];
    double V_2 = exp(-a*x_2 - b*tfinal)*Result[M][i + 1];
    double Vapprox = ((S_2 - data.S0)*V_1 + (data.S0 - S_1)*V_2) / (S_2 - S_1);
    
    return Vapprox;
}
