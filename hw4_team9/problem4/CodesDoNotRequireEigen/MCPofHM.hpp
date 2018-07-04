//
//  MCPofHM.hpp
//  Heston Model: Monte Carlo Pricing Heston Model
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef MCPofHM_hpp
#define MCPofHM_hpp

#include <iostream>
#include <tuple>
#include "OptionData.hpp"

typedef unsigned long int ULI;
typedef long double LDBL;

class HM
{ // Heston Model
    
private:
    option data;   // Encapulated data
  
public:
    // (0) Specification of option data
    void init();   // Initialization with arbitrary values
    void SetOption(double S0, double K, double T, double v, double q, double r);
    
    // (1) Constructor and destructor
    HM();          // Default constructor
    virtual ~HM(); // Destructor
    
    // (2) Monte Carlo pricing of Basket option
    double MC_HM(ULI n, ULI m, char PutCall, double lambda, double vbar, double eta, double rho);
    double implied_vol(double V, char PutCall);
};

#endif /* MCPofHM_hpp */
