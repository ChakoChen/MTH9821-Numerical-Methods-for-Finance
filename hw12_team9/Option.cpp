#include "Option.hpp"

Option::Option() {};

Option::Option(const Option & other) : S0(other.S0), T(other.T), K(other.K), r(other.r), sigma(other.sigma), q(other.q) {};

Option::Option(double S0, double T, double K, double r, double sigma, double q): S0(S0), T(T), K(K), r(r), sigma(sigma), q(q) {};

Option::~Option() {};

Option& Option::operator = (const Option & other) {

	if (this == &other) {
		return *this;
	}

	r = other.r, T = other.T, K = other.K;

	S0 = other.S0, sigma = other.sigma, q = other.q;

	return *this;
}

//get initial value S0
double Option::SpotPrice() {
	return S0;
}

//get initial value of volatility
double Option::Volatility() {
	return sigma;
}

//get time horizon

double Option::Maturity() {
	return T;
}
//get strike price
double Option::Strike() {
	return K;
}
//get dividend
double Option::Dividend() {
	return q;
}
//get interest rate
double Option::Interest() {
	return r;
}
