//
//  OptionPricing_FD.hpp
//  AmericanOption_FD
//
//  Finite difference methods for option pricing
//
//  Created by Changheng Chen on 10/15/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include <vector>
#include <tuple>
#include "OptionData.hpp"
using namespace std;

#ifndef OptionPricing_FD_hpp
#define OptionPricing_FD_hpp

#include <iostream>

class FDM_OP
{ // Finite difference methods for option pricing
    
private:
    option data; // Encapulated data
    
public:
    // (0) Specification of option data
    void init(); // Initialization with arbitrary values
    void SetOption(double S0, double K, double T, double v, double q, double r, char AmerEuro, char PutCall);
    
    // (1) Constructor and destructor
    FDM_OP();          // Default constructor
    virtual ~FDM_OP(); // Destructor
    
    // (2) Finite difference method
    tuple<vector<vector<double>>, vector<double>, vector<double>> FDM(double alpha, double omega, int M);
    double FDM(double alpha, int M);
};

#endif /* OptionPricing_FD_hpp */
