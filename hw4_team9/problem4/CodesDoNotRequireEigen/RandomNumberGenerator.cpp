//
//  RandomNumberGenerator.cpp
//  RandomNumberGenerator
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "RandomNumberGenerator.hpp"
#include <cmath>

// Constructor and destructor
RNG::RNG() {};
RNG::~RNG(){};

// Linear Congruential Generator
std::vector<LDBL> RNG::LCG(ULI size)
{
    std::vector<LDBL> rand(size);
    
    ULI a = 39373;
    ULI c = 0;
    ULI k = pow(2,31)-1;
    ULI x0 = 1;
    ULI x;

    x = x0;
    for (ULI i=0; i<size; i++)
    {
        x = (a*x+c) % k;
        rand[i] = (LDBL)x/(LDBL)k;
    }
    
    return rand;
}

// Inverse Transform Method
std::vector<LDBL> RNG::ITM(ULI size)
{
    std::vector<LDBL> rand(size);
    std::vector<LDBL> u = LCG(size);
    LDBL x, y, r;
    
    LDBL a0 = 2.50662823884;
    LDBL a1 = -18.61500062529;
    LDBL a2 = 41.39119773534;
    LDBL a3 = -25.44106049637;
    LDBL b0 = -8.47351093090;
    LDBL b1 = 23.08336743743;
    LDBL b2 = -21.06224101826;
    LDBL b3 = 3.13082909833;
    LDBL c0 = 0.3374754822726147;
    LDBL c1 = 0.9761690190917186;
    LDBL c2 = 0.16079797149118209;
    LDBL c3 = 0.0276438810333863;
    LDBL c4 = 0.0038405729373609;
    LDBL c5 = 0.0003951896511919;
    LDBL c6 = 0.0000321767881768;
    LDBL c7 = 0.0000002888167364;
    LDBL c8 = 0.0000003960315187;
    
    for (ULI i = 0; i<size; i++)
    {
        y = u[i]-0.5;
        if (std::abs(y)<0.42)
        {
            r = y*y;
            x = y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1);
        }
        else
        {
            r = u[i];
            if (y>0) r = 1-u[i];
            r = log(-log(r));
            x = c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
            if (y<0) x = -x;
        }
        rand[i] = x;
    }
    
    return rand;
}

// Acceptance-Rejection Method
std::vector<LDBL> RNG::ARM(ULI size)
{ // (# of LCG)/(# of ARM) > 3.95 (i.e. >395 LCG numbers can generate 100 BMM numbers)
    
    std::vector<LDBL> rand(size);
    std::vector<LDBL> u = LCG(size*4);
    LDBL X;
   
    ULI n = 0;
    for (ULI i=0; i<u.size()-2; i+=3)
    {
        X = -log(u[i]);
        if (u[i+1]>exp(-0.5*(X-1)*(X-1)))
            continue;
        else
        {
            if (u[i+2]<=0.5) X = -X; 
        }
        
        rand[n] = X;
        n++;
        
        if (n>=size) break;
    }
    
    return rand;
}

// Box-Muller Method
std::vector<LDBL> RNG::BMM(ULI size)
{ // (# of LCG)/(# of BMM) ~= 1.27 (i.e. ~127 LCG numbers can generate 100 BMM numbers)
    
    std::vector<LDBL> rand(size);
    std::vector<LDBL> u = LCG(size*1.5);
    LDBL u1, u2, X, Y, Z1, Z2;
    
    ULI n = 0;
    for (ULI i=0; i<u.size()-1; i+=2)
    {
        u1 = 2*u[i]-1; u2 = 2*u[i+1]-1;
        X = u1*u1 + u2*u2;
        if (X>1) continue;
        
        Y = sqrt(-2*log(X)/X);
        Z1 = u1*Y; Z2 = u2*Y;
        
        rand[n] = Z1; rand[n+1] = Z2;
        n+=2;
        
        if (n>=size) break;
    }
    
    return rand;
}
