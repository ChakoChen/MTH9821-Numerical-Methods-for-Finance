//BarrierOption
#ifndef BarrierOption_HPP
#define BarrierOption_HPP

#include <iostream>
#include <vector>
#include <string>
using namespace std;

//a struct to combine the result of TriDiagonal LU decomposition
struct LUResult {
	vector<vector<double>> Lmatrix;
	vector<vector<double>> Umatrix;
};
//To do Tridiagonal LU decomposition on A;
LUResult Tridiagonal_LU(vector<vector<double>> A);

//pdf of standard normal
double pdf(double x);

//cdf of standard normal
double cdf(double x);




class BarrierOption {
private:
	double S;
	double K;
	double B;
	double T;
	double q;
	double r;
	double sigma;
public:
	BarrierOption(double s, double k, double b, double t, double qq, double rr, double sig)
		:S(s), K(k), B(b), T(t), q(qq), r(rr), sigma(sig) {};
	~BarrierOption() {};

	//BlackScholes for Call
	vector<double> Black_Scholes_Vanilla(string type);
	//BlackScholes for Barrier Option
	vector<double> Black_Scholes_Barrier();

	//Domain Discretization
	vector<double> DomainDiscret(double alpha, int M);

	//ForwardEuler
	//if Umatrix = ¡°Yes¡±: Print Umatrix
	vector<double> ForwardEuler(double alpha, int M, string Umatrix);

	//BackwardEuler_LU
	//if Umatrix = ¡°Yes¡±: Print Umatrix
	vector<double> BackwardEuler_LU(double alpha, int M, string Umatrix);

	//CrankNicloson_SOR
	vector<double> CrankNicolson_SOR(double alpha, int M, double w, double tol);
};

#endif