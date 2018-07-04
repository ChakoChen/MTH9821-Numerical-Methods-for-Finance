//
//  MCPofBO.hpp
//  BasketOption: Monte Carlo Pricing of Basket Option
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef MCPofBO_hpp
#define MCPofBO_hpp

#include <iostream>
#include <tuple>
#include "OptionData.hpp"

typedef unsigned long int ULI;
typedef long double LDBL;

class BO
{ // Basket Option class
    
private:
    basket data;   // Encapulated data
  
public:
    // (0) Specification of option data
    void initBasket();   // Initialization with arbitrary values
    void SetBasket(double S01, double S02, double v1, double v2, double rho, double K, double T, double q, double r);
    
    // (1) Constructor and destructor
    BO();          // Default constructor
    virtual ~BO(); // Destructor
    
    // (2) Monte Carlo pricing of Basket option
    double MC_BO(ULI N, char PutCall);
    double MC_Path_BO(ULI n, ULI m, char PutCall);
};

#endif /* MCPofBO_hpp */
