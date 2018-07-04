//
//  BinomialTreePricer.hpp
//  BinomialTreeII
//
//  Created by Changheng Chen on 9/3/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef BinomialTreePricer_hpp
#define BinomialTreePricer_hpp

#include <iostream>
#include <vector>
#include <array>

// (1) Binomial tree method
double BT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double BTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

// (2) Averaged binomial tree method
double ABT(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double ABTGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

// (3) BBS
double BBS(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double BBSGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

// (4) BBSR
double BBSR(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer);
double BBSRGreeks(int N, double S0, double K, double T, double q, double r, double v, char PutCall, char EuroAmer, char Greek);

#endif /* BinomialTreePricer_hpp */
