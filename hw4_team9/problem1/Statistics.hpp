//
//  Statistics.hpp
//  VarianceReductionTechniques
//
//  Functions for calcuculating the mean of a vector and covariance of two vectors
//
//  Created by Changheng Chen on 9/20/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#ifndef Statistics_hpp
#define Statistics_hpp

#include <iostream>
#include <vector>
#include <numeric>

double mean(std::vector<double> X)
{
    double Xbar = std::accumulate(X.begin(), X.end(), 0.0);
    return Xbar/X.size();
}

double covariance(std::vector<double> X, std::vector<double> Y)
{
    double Xbar = mean(X);
    double Ybar = mean(Y);
    double cov  = .0;
    
    for (ULI i=0; i<X.size(); i++)
    {
        cov += (X[i]-Xbar)*(Y[i]-Ybar);
    }
    
    return cov;
}

#endif /* Statistics_hpp */
