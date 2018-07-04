//
//  main.cpp
//  AmericanOption_FD
//
//  Created by Changheng Chen on 10/14/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "OptionPricing_FD.hpp"
#include "BlackScholes.hpp"
#include <iomanip>
#include <iostream>

using namespace std;

void P1()
{ // Pricing European Put Option
    
    double S0 = 42;
    double K  = 40;
    double T  = 0.75;
    double v  = 0.32;
    double q  = 0.02;
    double r  = 0.04;
    char AmerEuro = 'A';
    char PutCall = 'P';
    
    // Get solution u for alpha=0.45 and M=4
    int M = 4; double alpha = 0.45; double omega = 1.2;
    FDM_OP option;
    option.SetOption(S0, K, T, v, q, r, AmerEuro, PutCall);
    tuple<vector<vector<double>>, vector<double>, vector<double>> output = option.FDM(alpha, omega, M);
    vector<vector<double>> u = get<0>(output);
    
    // Print u for M=4 and N=11
    cout << setprecision(12);
    cout << "m x_left x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10 x_right" << endl;
    for (int m=0; m<u[0].size(); m++)
    {
        cout << m << ", ";
        for (int n=0; n<u.size(); n++)
        {
            cout << u[n][m] << ", ";
        }
        cout << endl;
    }
    
    // Change M and get error statistics and the Greeks
    cout << "\nM error_pointwise Ratio(error_pointwise) error_pointwise_2 Ratio(error_pointwise2) Delta Gamma Theta VarRed VarRed(error_pointwise)" << endl;
    vector<double> statistics_old, statistics_new, greeks;
    for (int M=4; M<=256; M*=4)
    {
        output = option.FDM(alpha, omega, M);
        statistics_new = get<1>(output);
        greeks = get<2>(output);
        
        cout << M << ", ";
        if (M==4)
        {
            cout << statistics_new[0] << ", , " << statistics_new[1] << ", , "
                 << greeks[0] << ", " << greeks[1] << ", " << greeks[2] << ", "
                 << statistics_new[2] << ", " << statistics_new[3] << endl;
        }
        else
        {
            cout << statistics_new[0] << ", " << statistics_new[0]/statistics_old[0] << ", "
                 << statistics_new[1] << ", " << statistics_new[1]/statistics_old[1] << ", "
                 << greeks[0] << ", " << greeks[1] << ", " << greeks[2] << ", "
                 << statistics_new[2] << ", " << statistics_new[3] << endl;
        }
        statistics_old = statistics_new;
    }
    
    alpha = 5.;
    // Change M and get error statistics and the Greeks
    cout << "\nM error_pointwise Ratio(error_pointwise) error_pointwise_2 Ratio(error_pointwise2) Delta Gamma Theta VarRed VarRed(error_pointwise)" << endl;
    for (int M=4; M<=256; M*=4)
    {
        output = option.FDM(alpha, omega, M);
        statistics_new = get<1>(output);
        greeks = get<2>(output);
        
        cout << M << ", ";
        if (M==4)
        {
            cout << statistics_new[0] << ", , " << statistics_new[1] << ", , "
            << greeks[0] << ", " << greeks[1] << ", " << greeks[2] << ", "
            << statistics_new[2] << ", " << statistics_new[3] << endl;
        }
        else
        {
            cout << statistics_new[0] << ", " << statistics_new[0]/statistics_old[0] << ", "
            << statistics_new[1] << ", " << statistics_new[1]/statistics_old[1] << ", "
            << greeks[0] << ", " << greeks[1] << ", " << greeks[2] << ", "
            << statistics_new[2] << ", " << statistics_new[3] << endl;
        }
        statistics_old = statistics_new;
    }
}


int main()
{
    P1();
    
    return 0;
}
