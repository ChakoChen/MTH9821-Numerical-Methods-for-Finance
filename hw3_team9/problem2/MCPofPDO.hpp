//
//  MCPofPDO.hpp
//  PathDependentOption: Monte Carlo Pricing of a Path-Dependent Option
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef MCPofPDO_hpp
#define MCPofPDO_hpp

#include <iostream>
#include "OptionData.hpp"

typedef unsigned long int ULI;
typedef long double LDBL;

class PDO
{ // Path-Dependent Option class
    
private:
    option data; // Encapulated data
  
public:
    // (0) Specification of option data
    void init();                  // Initialization with arbitrary values
    void SetOption(double S0, double K, double B, double T, double v, double q, double r);
    
    // (1) Constructors, destructor, and assignment operator
    PDO();                        // Default constructor
    virtual ~PDO();               // Destructor
    
    // (2) Pricing of option
    double exact_Cdao();          // Closed formula
    double MC_Cdao(ULI n, ULI m); // Monte Carlo pricing
};

#endif /* MCPofPDO_hpp */
