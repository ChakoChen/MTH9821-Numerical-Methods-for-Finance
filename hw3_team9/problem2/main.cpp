//
//  main.cpp
//  PathDependentOption
//
//  Created by Changheng Chen on 9/16/17.
//  Copyright © 2017 Changheng Chen. All rights reserved.
//

#include "MCPofPDO.hpp"
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <iomanip>

double S0 = 39;
double K  = 39;
double B  = 35;
double T  = 0.75;
double v  = 0.25;
double q  = 0.01;
double r  = 0.02;

int main()
{
    ULI m, n, Nk;
    
    PDO MyOption;
    MyOption.SetOption(S0, K, B, T, v, q, r);
    double Cdao = MyOption.exact_Cdao();
    double MC_C;
    
    //============ 3.1 Fixed m ============
    std::cout << "m, n, V(n), |Cdao-V(n)|" << std::endl;
    m = 200;
    for (int k=0; k<10; k++)
    {
        Nk = 10000*pow(2,k);
        MC_C = MyOption.MC_Cdao(Nk/m, m);
        
        std::cout << std::setprecision(9);
        std::cout << m << ", " << Nk/m << ", " << MC_C << ", "
                  << std::abs(Cdao-MC_C) << std::endl;
    }
    
    //=========  3.2 Optimal m, n =========
    std::cout << "\nmk, nk, V(n), |Cdao-V(n)|" << std::endl;
    for (int k=0; k<10; k++)
    {
        Nk = 10000*pow(2,k);
        m = (ULI)ceil(pow((double)Nk,1./3.)*pow(T,2./3.));
        n = (ULI)floor(Nk/m);
        MC_C = MyOption.MC_Cdao(n, m);
        
        std::cout << m << ", " << n << ", " << MC_C << ", "
        << std::abs(Cdao-MC_C) << std::endl;
    }

    return 0;
}
