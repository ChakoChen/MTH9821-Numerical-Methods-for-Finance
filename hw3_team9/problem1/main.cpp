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

double S0 = 41;
double K  = 42;
double T  = 0.75;
double v  = 0.20;
double q  = 0.01;
double r  = 0.03;
char PutCall;

int main()
{
    
    cout << "Please enter: 'C' for European call, 'P' for European put: ";
    cin >> PutCall;
    cout << setprecision(12);
    
    std::tuple<double, double, double> OptionGreeks;
    double V_BS, delta_BS, vega_BS;
    ULI N;
    EO MyOption;
    MyOption.SetOption(S0, K, T, v, q, r);

    // (1) Black-Scholes option ad Greeks
    V_BS = BlackScholes(S0, K, T, v, q, r, PutCall);
    delta_BS = BlackScholesGreeks(S0, K, T, v, q, r, PutCall,'D');
    vega_BS = BlackScholesGreeks(S0, K, T, v, q, r, PutCall,'V');

    cout << "Black-Scholes: V, delta, vega \n";
    cout << V_BS << ", " << delta_BS << ", " << vega_BS << "\n\n";
    
    // (2) Monte Carlo pricing and Greeks estimations
    cout << "Monte Carlo: N, V, sqrt(N)|V_BS-V|, delta, sqrt(N)|delta_BS-delat|, vega, sqrt(N)|vega_BS-vega| \n";
    for (int k=0; k<10; k++)
    {
        N = 10000*pow(2,k);
        OptionGreeks = MyOption.MC_EO(N, PutCall);
        cout << N << ", " << get<0>(OptionGreeks) << ", " << sqrt(N)*abs(V_BS-get<0>(OptionGreeks)) << ", " << get<1>(OptionGreeks) << ", " << sqrt(N)*abs(delta_BS-get<1>(OptionGreeks)) << ", " << get<2>(OptionGreeks) << ", " << sqrt(N)*abs(vega_BS-get<2>(OptionGreeks)) << "\n";
    }
       
    return 0;
}
