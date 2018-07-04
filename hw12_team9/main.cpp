//
// Longstaff_Schwartz.cpp and compare results with binomial tree method
//
#include "PathProcess.hpp"
#include <iostream>
#include <iomanip>
#include <fstream> 

using namespace std;

double S = 40;
double K = 45;
double T = 0.5;
double r = 0.05;
double v = 0.23;
double q = 0.0;

int main()
{
	shared_ptr<AmericanPut> AmerOptr = shared_ptr<AmericanPut>(new AmericanPut(S, T, K, r, v, q));

	int M = 1000, N = 100000, B = 10;
    double res;
	
    cout << "Power Basis:\n";
	for (int i=3; i<=B; i++)
    {
        res = PathProcess(AmerOptr, M, N, i, "Box Muller", "Power Basis");
        cout << res << "\n";
	}

    cout << "Weighted Laguerre Polynomials Basis:\n";
	for (int i=3; i<=B; i++)
    {
        res = PathProcess(AmerOptr, M, N, i, "Box Muller", "Laguerre Polynomials");
        cout << res << "\n";
	}
    
	return 0;
}
