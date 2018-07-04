#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <algorithm>
#include "OptionPricing.hpp"
using namespace std;

//pdf of standard normal
double pdf(double x) {
	return exp(-pow(x, 2) / 2) / (sqrt(2 * M_PI));
}

//cdf of standard normal
double cdf(double x) {
	return  0.5 + erf(x / sqrt(2)) / 2;
}

//BlackScholes Value and Greeks(Delta, Gamma, Theta, Vega)
vector<double> OptionPricing::BlackScholes(string type) {
	double d1 = (log(S / K) + (r - q + pow(sigma, 2)*0.5)*T) / (sigma*sqrt(T));
	double d2 = d1 - sigma*sqrt(T);

	double calltheta = q*S*exp(-q*T)*cdf(d1) - r*K*exp(-r*T)*cdf(d2) -
		K*exp(-r*T)*pdf(d2)*sigma / (2 * sqrt(T));

	vector<double> Result;
	if (type.compare("PUT") == 0) {
		double P = K*exp(-r*T)*cdf(-d2) - S*exp(-q*T)*cdf(-d1);			//PValue
		double DeltaP = exp(-q*T)*(cdf(d1) - 1);							//PDelta
		double GammaP = exp(-q*T)*pdf(d1) / (S*sigma*sqrt(T));			//PGamma
		double ThetaP = calltheta - q*S*exp(-q*T) + r*K*exp(-r*T);		//PTheta
		double VegaP = S*exp(-q*T)*pdf(d1)*sqrt(T);
		Result = { P, DeltaP, GammaP, ThetaP, VegaP };
	}
	else if (type.compare("CALL") == 0) {
		double C = S*exp(-q*T)*cdf(d1) - K*exp(-r*T)*cdf(d2);			//Cvalue
		double DeltaC = exp(-q*T)*cdf(d1);								//CDelta
		double GammaC = exp(-q*T)*pdf(d1) / (S*sigma*sqrt(T));			//CGamma
		double ThetaC = calltheta;										//CTheta
		double VegaC = S*exp(-q*T)*pdf(d1)*sqrt(T);
		Result = { C, DeltaC, GammaC, ThetaC, VegaC };
	}
	return Result;
}

