//
//  MCPofEO.hpp
//  PathDependentOption: Monte Carlo Pricing of European Option
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef MCPofEO_hpp
#define MCPofEO_hpp

#include <iostream>
#include <tuple>
#include "OptionData.hpp"

typedef unsigned long int ULI;
typedef long double LDBL;

class EO
{ // Path-Dependent Option class
    
private:
    option data; // Encapulated data
  
public:
    // (0) Specification of option data
    void init();   // Initialization with arbitrary values
    void SetOption(double S0, double K, double T, double v, double q, double r);
    
    // (1) Constructors, destructor, and assignment operator
    EO();          // Default constructor
    virtual ~EO(); // Destructor
    
    // (2) Pricing of option
    std::tuple<double, double, double> MC_EO(ULI N, char PutCall); // Monte Carlo pricing of European option
};

#endif /* MCPofEO_hpp */
