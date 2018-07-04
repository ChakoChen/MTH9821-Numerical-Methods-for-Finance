#ifndef OPTION_HPP
#define OPTION_HPP

#include <iostream>
#include <vector>
#include <string>
using namespace std;

//a struct to combine the result of TriDiagonal LU decomposition
struct LUResult {
	vector<vector<double>> Lmatrix;
	vector<vector<double>> Umatrix;
};
//a struct to combine the result of the Pricing Function
struct PriceResult {
	vector<vector<double>> Umatrix;
	vector<double> PointwiseConvergence;
};
//To do Tridiagonal LU decomposition on A;
LUResult Tridiagonal_LU(vector<vector<double>> A);
//cdf of standard normal
double cdf(double x);
//pdf of standard normal
double pdf(double x);


class Option {
private:
	double S;
	double K;
	double T;
	double r;
	double q;
	double sigma;
	string type;		//"PUT" or "CALL"
public:
	//constructor
	Option(double S, double K, double T, double r, double q, double sigma, string type):
		S(S), K(K), T(T), r(r), q(q), sigma(sigma), type(type){};
	//destructor
	~Option() {};

	//member functions
	//return be values calculated by Black_Scholes
	vector<double> BlackScholes();
	//return be values calculated by Backward Euler Method
	PriceResult BackWardEuler(int M, double alpha_temp);
	//return be values calculated by CrankNicolson Method
	PriceResult CrankNicolson(int M, double alpha_temp);
};



#endif