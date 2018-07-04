//
//  main.cpp
//  TrinomialTreeIII
//
//  Created by Changheng Chen on 9/10/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "TrinomialTreePricer.hpp"
#include "BlackScholes.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <stdio.h>
#include <fstream> // file
#include <iomanip> // setprecision

using namespace std;

// Trinomial tree parameters
int N = 100;          // Initial value of the number of time intervals
double S0 = 41;       // Spot price
double K = 39;        // Strike price
double T = 1;         // Maturity
double q = 0.005;     // Dividends
double r =  0.03;     // Risk free interest rates
double v =  0.25;     // Volatility
char PutCall = 'P';   // 'C' for call, 'P' for put
char EuroAmer = 'E';  // 'A' for American option, 'E' for European option
int method = 4;       // 1 = Trinomial tree; 2 = Averaged Trinomial Tree; 3 = 'TBS'; 4 ='TBSR'
std::string filename; // Output filename

int main()
{
    // Specify the methods to be used and the output file name 
    cout << "Please enter a file name to write: ";
    std::getline(std::cin, filename);
    cout << "Please enter: 'A' for American option, 'E' for European option ";
    cin >> EuroAmer;
    cout << "Please enter a number for: 1 = Trinomial tree; 2 = Averaged Trinomial Tree; 3 = 'TBS'; 4 ='TBSR' ";
    cin >> method;

    // Black-Scholes result
    double V_BS = BlackScholes(S0, K, T, v, q, r, PutCall);
    double Delta_BS = PutDelta(T, K, v, r, q, S0);
    double Gamma_BS = Gamma(T, K, v, r, q, S0);
    double Theta_BS = PutTheta(T, K, v, r, q, S0);
    
    cout << std::setprecision(8) << "V_BS:      " << V_BS << endl;     //  V_BS =  2.6112396
    cout << std::setprecision(8) << "Delta_BS: "  << Delta_BS << endl; // Delta = -0.33373031
    cout << std::setprecision(8) << "Gamma_BS:  " << Gamma_BS << endl; // Gamma =  0.035382198
    cout << std::setprecision(8) << "Theta_BS: "  << Theta_BS << endl; // Theta = -1.4382604
    
    // Compute the option using the 4 methods
    vector<vector<double>> V(10); //Vectors to save the 10 required variables
    double price, delta, gamma, theta;
    
    for (N=10; N<=1280; N=N*2)
    {
        cout << "N = " << N << endl;
        
        price = Pricer(N, S0, K, T, q, r, v, PutCall, EuroAmer, method);
        delta = Greeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'D', method);
        gamma = Greeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'G', method);
        theta = Greeks(N, S0, K, T, q, r, v, PutCall, EuroAmer,'T', method);
        
        V[0].push_back(price);
        V[1].push_back(abs(price-V_BS));
        V[2].push_back(N*abs(price-V_BS));
        V[3].push_back(N*N*abs(price-V_BS));
        V[4].push_back(delta);
        V[5].push_back(abs(delta-Delta_BS));
        V[6].push_back(gamma);
        V[7].push_back(abs(gamma-Gamma_BS));
        V[8].push_back(theta);
        V[9].push_back(abs(theta-Theta_BS));
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
