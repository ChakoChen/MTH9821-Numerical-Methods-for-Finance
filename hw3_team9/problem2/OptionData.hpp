//
//  OptionData.hpp
//  PathDependentOption
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef OptionData_hpp
#define OptionData_hpp

#include <iostream>

struct option
{
    double S0; // Spot price
    double K;  // Strike price
    double B;  // Barrier
    double T;  // Maturity
    double v;  // volatility
    double q;  // dividends
    double r;  // risk-free risk
};

#endif /* OptionData_hpp */
