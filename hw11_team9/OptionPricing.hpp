#ifndef OPTION_PRICING
#define OPTION_PRICING

#include <iostream>
#include <vector>
#include <string>
using namespace std;

class OptionPricing {
private:
	double S;
	double K;
	double T;
	double r;
	double q;
	double sigma;

public:
	OptionPricing() {};
	OptionPricing(double s, double k, double t, double rr, double qq, double ss) :
		S(s), K(k), T(t), r(rr), q(qq), sigma(ss) {};
	virtual ~OptionPricing() {};

	double getS() { return S; };
	double getK() { return K; };
	double getT() { return T; };
	double getr() { return r; };
	double getq() { return q; };
	double getsigma() { return sigma; };

	//BlackScholes Value and Greeks
	vector<double> BlackScholes(string type);

};

//pdf of standard normal
double pdf(double x);

//cdf of standard normal
double cdf(double x);

#endif // !OPTION_PRICING