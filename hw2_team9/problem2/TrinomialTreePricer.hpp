//
//  TrinomialTreePricer.hpp
//  TrinomialTreeII
//
//  Created by Changheng Chen on 9/10/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef TrinomialTreePricer_hpp
#define TrinomialTreePricer_hpp

#include "BlackScholes.hpp"

// (0) Define functions to call the four methods with a switch key "method"
double Pricer(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, int method);
double Greeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek, int method);

// (1) Trinomial tree method
double TT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double TTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

// (2) Averaged Trinomial tree method
double ATT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double ATTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

// (3) TBS
double TBS(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double TBSGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

// (4) TBSR
double TBSR(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double TBSRGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

#endif /* TrinomialTreePricer_hpp */
