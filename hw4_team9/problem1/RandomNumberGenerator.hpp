//
//  RandomNumberGenerator.hpp
//  RandomNumberGenerator
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef RandomNumberGenerator_hpp
#define RandomNumberGenerator_hpp

#include <iostream>
#include <vector>

typedef unsigned long int ULI;
typedef long double LDBL;

class RNG
{
public:
    RNG();                            // Default constructor
    virtual ~RNG();                   // Destructor

    std::vector<LDBL> LCG(ULI size);  // Linear Congruential Generator
    std::vector<LDBL> ITM(ULI size);  // Inverse Transform Method
    std::vector<LDBL> ARM(ULI size);  // Acceptance-Rejection Method
    std::vector<LDBL> BMM(ULI size);  // Box-Muller Method
};

#endif /* RandomNumberGenerator_hpp */
