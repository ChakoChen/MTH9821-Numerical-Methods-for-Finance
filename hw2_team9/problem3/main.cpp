//
//  main.cpp
//  TrinomialTreeIII
//
//  Created by Changheng Chen on 9/10/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "BlackScholes.hpp"
#include "BinomialTreePricer.hpp"
#include "TrinomialTreePricer.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <stdio.h>
#include <iostream> 
#include <fstream>  // file
#include <iomanip>  // setprecision

using namespace std;

// Trinomial tree parameters
int N = 100;           // Initial value of the number of time intervals
double S0 = 41;        // Spot price
double K = 39;         // Strike price
double T = 1;          // Maturity
double q = 0.005;      // Dividends
double r =  0.03;      // Risk free interest rates
double v =  0.25;      // Volatility
char PutCall = 'P';    // 'C' for call, 'P' for put
char EuroAmer = 'A';   // 'A' for American option, 'E' for European option
char VarRed = 'Y';     // 'Y' for variance reduction, 'N' for without variance reduction
int method = 4;        // 1 = Trinomial tree; 2 = Averaged Trinomial Tree; 3 = 'TBS'; 4 ='TBSR'
std::string filename;  // Output filename

int main()
{
    // Specify the methods to be used and the output file name 
    cout << "Please enter a file name to write: ";
    std::getline(std::cin, filename);
    cout << "Please enter: 'Y' for variance reduction, 'N' for without variance reduction ";
    cin >> VarRed;
    cout << "Please enter a number for: 1 = Trinomial tree; 2 = Averaged Trinomial Tree; 3 = 'TBS'; 4 ='TBSR' ";
    cin >> method;

    /*// Run ATT only one time to approximate exact solutions
    N = 10000;
    double V_exact = ABT(N, S0, K, T, q, r, v, PutCall, EuroAmer);
    double Delta_exact = ABTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'D');
    double Gamma_exact = ABTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'G');
    double Theta_exact = ABTGreeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'T');
    
    cout << std::setprecision(16) << "V_exact:      " << V_exact << endl;
    cout << std::setprecision(16) << "Delta_exact: "  << Delta_exact << endl;
    cout << std::setprecision(16) << "Gamma_exact:  " << Gamma_exact << endl;
    cout << std::setprecision(16) << "Theta_exact: "  << Theta_exact << endl;
    */
    
    double V_exact =      2.67707291925769;
    double Delta_exact = -0.3453953479455988;
    double Gamma_exact =  0.03738936277481288;
    double Theta_exact = -1.529771437526925;

    // Compute the option using the 4 methods
    vector<vector<double> > V(10); //Vectors to save the 10 required variables
    double price, delta, gamma, theta;
    
    for (N=10; N<=1280; N=N*2)
    {
        cout << "N = " << N << endl;
        
        switch(VarRed)
        {
            case('N'):
            {
                price = Pricer(N, S0, K, T, q, r, v, PutCall, EuroAmer, method);
                delta = Greeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'D', method);
                gamma = Greeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'G', method);
                theta = Greeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'T', method);
                break;
            }
            case('Y'):
            {
                price = PricerRD(N, S0, K, T, q, r, v, PutCall, EuroAmer, method);
                delta = GreeksRD(N, S0, K, T, q, r, v, PutCall, EuroAmer,'D', method);
                gamma = GreeksRD(N, S0, K, T, q, r, v, PutCall, EuroAmer,'G', method);
                theta = GreeksRD(N, S0, K, T, q, r, v, PutCall, EuroAmer,'T', method);
                break;
            }
        }
        
        V[0].push_back(price);
        V[1].push_back(abs(price-V_exact));
        V[2].push_back(N*abs(price-V_exact));
        V[3].push_back(N*N*abs(price-V_exact));
        V[4].push_back(delta);
        V[5].push_back(abs(delta-Delta_exact));
        V[6].push_back(gamma);
        V[7].push_back(abs(gamma-Gamma_exact));
        V[8].push_back(theta);
        V[9].push_back(abs(theta-Theta_exact));
    }
    
    // Change the file name for differenc method!
    ofstream f(filename.c_str(), std::ios::out | std::ios::trunc);
    ostream_iterator<double> output_iterator(f, " ");
    f.precision(8);
    f.setf(ios::fixed);
    for (int i=0; i<10; i++)
    {
        copy(V[i].begin(), V[i].end(), output_iterator);
        f << '\n';
    }
    f.close();
    
    return 0;
}
