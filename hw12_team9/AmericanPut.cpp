#include "AmericanPut.hpp"

AmericanPut::AmericanPut() {};

AmericanPut::AmericanPut(const AmericanPut & other) : Option(other){};

AmericanPut::AmericanPut(double S0, double T, double K, double r, double sigma, double q) : Option(S0, T, K, r, sigma, q) {};

AmericanPut::AmericanPut(double S0, double T, double K, double r, double sigma, double q, double exact) : Option(S0, T, K, r, sigma, q), V_exact(exact) {};

AmericanPut::~AmericanPut(){};

AmericanPut& AmericanPut::operator = (const AmericanPut & option)
{

	if (this==&option)
        return *this;
	Option::operator=(option);
	return *this;
}

std::string AmericanPut::OptionType()
{
	return "American Put";
}

double AmericanPut::Price()
{
	return V_exact;
}
