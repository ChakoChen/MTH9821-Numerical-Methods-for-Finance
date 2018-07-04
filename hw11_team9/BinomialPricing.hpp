#ifndef BINOMIAL_PRICING
#define BINOMIAL_PRICING

#include "OptionPricing.hpp"
using namespace std;

class Binomial : public OptionPricing {
public:
	Binomial() {};
	Binomial(double s, double k, double t, double rr, double qq, double ss) :OptionPricing(s, k, t, rr, qq, ss) {};
	~Binomial() {};

	vector<double> American_NoDiv(string type, int N);
	vector<double> European_Nodiv(string type, int N);

	vector<double> European_PropDiv(string type, int N);
	vector<double> American_PropDiv(string type, int N, double t_div);
	vector<double> European_PropDiv_Mult(string type, int N, int count);
	vector<double> American_PropDiv_Mult(string type, int N, vector<double> times);

	vector<double> European_FixDiv(string type, int N, double t_div, double v_div);
	vector<double> American_FixDiv(string type, int N, double t_div, double v_div);
	vector<double> European_FixDiv_Mult(string type, int N, vector<double> times, double v_div);
	vector<double> American_FixDiv_Mult(string type, int N, vector<double> times, double v_div);


	vector<double> European_Complex(string type, int N, vector<double> time, vector<double> value);
	vector<double> American_Complex(string type, int N, vector<double> time, vector<double> value);
};



#endif // !BINOMIAL_PRICING
