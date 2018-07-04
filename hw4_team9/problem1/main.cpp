//
//  main.cpp
//  PathDependentOption
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <cmath>
#include <iomanip>
#include <iostream>
#include "MCPofEO.hpp"
#include "BlackScholes.hpp"
using namespace std;

double S0 = 50;
double K  = 55;
double T  = 0.75;
double v  = 0.30;
double q  = 0.00;
double r  = 0.04;
char PutCall = 'P';

int main()
{
    
    // cout << "Please enter: 'C' for European call, 'P' for European put: ";
    // cin >> PutCall;
    cout << setprecision(12);

    double V_BS, V_MC; // Option price from Black-Scholes and Monte Carlo
    ULI N;             // Number of Monte Carlo simulations
    
    // (0) Creat option object
    EO MyOption; MyOption.SetOption(S0, K, T, v, q, r);

    // (1) Black-Scholes option
    V_BS = BlackScholes(S0, K, T, v, q, r, PutCall);
    cout << "Black-Scholes: V_BS = " << V_BS << "\n\n";
    
    // (2.1) Monte Carlo pricing with control variate technique
    cout << "Control Variate: N, V_CV, |V_BS-V_CV| \n";
    for (ULI k=0; k<10; k++)
    {
        N = 10000*pow(2,k);
        V_MC = MyOption.MC_EO_CV(N, PutCall);
        cout << N << ", " << V_MC << ", " <<  abs(V_BS-V_MC) << "\n";
    }
    cout << endl; 
    
    // (2.2) Monte Carlo pricing with antithetic variables
    cout << "Anthithetic Variables: N, V_AV, |V_BS-V_AV| \n";
    for (ULI k=0; k<10; k++)
    {
        N = 10000*pow(2,k);
        V_MC = MyOption.MC_EO_AV(N, PutCall);
        cout << N << ", " << V_MC << ", " <<  abs(V_BS-V_MC) << "\n";
    }
    cout << endl;
    
    // (2.3) Monte Carlo pricing with moment matching
    cout << "Moment Matching: N, V_MM, |V_BS-V_MM| \n";
    for (ULI k=0; k<10; k++)
    {
        N = 10000*pow(2,k);
        V_MC = MyOption.MC_EO_MM(N, PutCall);
        cout << N << ", " << V_MC << ", " <<  abs(V_BS-V_MC) << "\n";
    }
    cout << endl;
    
    // (2.4) Monte Carlo pricing with moment matching and contral variates
    cout << "Moment Matching and Control Variates: N, V_MM, |V_BS-V_MM| \n";
    for (ULI k=0; k<10; k++)
    {
        N = 10000*pow(2,k);
        V_MC = MyOption.MC_EO_MMandCV(N, PutCall);
        cout << N << ", " << V_MC << ", " <<  abs(V_BS-V_MC) << "\n";
    }
       
    return 0;
}
