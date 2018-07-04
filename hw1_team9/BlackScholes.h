//Black Scholes Pricer 

#ifndef BLACKSCHOLES_HPP
#define BLACKSCHOLES_HPP

#include <iostream>
#include <cmath>
#include <vector>

struct Option
{
	double K;
	double S0;
	double q;
	double vol;
	double r;
	double T;
};

double pi = 3.14159265358979323846;

//Calculate N(t)
double normalDistribution(double t)
{
		return 0.5*std::erfc(-t/std::sqrt(2));
}

std::vector<double> BSOption(const Option& op, const char & pc)
{
	double d1 = (std::log(op.S0 / op.K) + (op.r - op.q + op.vol*op.vol / 2.0)*(op.T)) / (op.vol*std::sqrt(op.T));
	double d2 = d1 - op.vol*std::sqrt(op.T);
	std::vector<double> result; //{delta,gamma,theta}
	double P;
	//price
	if(pc=='p')
	{//put option
		P = op.K*std::exp(-op.r*op.T)*normalDistribution(-d2) - op.S0*std::exp(-op.q*op.T)*normalDistribution(-d1);
	}
	else
	{//call option
		P = op.S0*std::exp(-op.q*op.T)*normalDistribution(d1)-op.K*std::exp(-op.r*op.T)*normalDistribution(d2);
	}
	result.push_back(P);
	//delta
	result.push_back( -std::exp(-op.q*op.T)*normalDistribution(-d1));
	//gamma
	result.push_back((std::exp(-op.q*op.T))/(op.S0*op.vol*std::sqrt(op.T))*(1)/(std::sqrt(2*pi))*std::exp(-d1*d1/2.0));
	//theta
	result.push_back(-(op.S0*op.vol*std::exp(-op.q*op.T)) / (2*std::sqrt(2*pi*op.T))*std::exp(-d1*d1/2.0) - op.q*op.S0*std::exp(-op.q*op.T)*normalDistribution(-d1) + op.r*op.K*std::exp(-op.r*op.T)*normalDistribution(-d2));
	
	return result;
}
#endif

