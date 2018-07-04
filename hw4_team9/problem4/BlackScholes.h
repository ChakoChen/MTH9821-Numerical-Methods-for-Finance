//Black Scholes Pricer 

#ifndef BLACKSCHOLES_HPP
#define BLACKSCHOLES_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "equationsolver.h"

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

//return {price, delta, gamma, theta}
std::vector<double> BSOption(const Option& op, const char & pc)
{
	double d1 = (std::log(op.S0 / op.K) + (op.r - op.q + op.vol*op.vol / 2.0)*(op.T)) / (op.vol*std::sqrt(op.T));
	double d2 = d1 - op.vol*std::sqrt(op.T);
	std::vector<double> result; //{delta,gamma,theta}
	if(pc=='p')
	{//put option
		//price
		result.push_back(op.K*std::exp(-op.r*op.T)*normalDistribution(-d2) - op.S0*std::exp(-op.q*op.T)*normalDistribution(-d1));
		//delta
		result.push_back(-std::exp(-op.q*op.T)*normalDistribution(-d1));
		//gamma
		result.push_back((std::exp(-op.q*op.T)) / (op.S0*op.vol*std::sqrt(op.T))*(1) / (std::sqrt(2 * pi))*std::exp(-d1*d1 / 2.0));
		//theta
		result.push_back(-(op.S0*op.vol*std::exp(-op.q*op.T)) / (2 * std::sqrt(2 * pi*op.T))*std::exp(-d1*d1 / 2.0) - op.q*op.S0*std::exp(-op.q*op.T)*normalDistribution(-d1) + op.r*op.K*std::exp(-op.r*op.T)*normalDistribution(-d2));
		//vega
		result.push_back(op.S0*std::exp(-op.q*op.T)*std::sqrt(op.T)*(1.0 / std::sqrt(2.0*pi))*std::exp(-d1*d1 / 2.0));
	}
	else
	{//call option
		//price
		result.push_back( op.S0*std::exp(-op.q*op.T)*normalDistribution(d1)-op.K*std::exp(-op.r*op.T)*normalDistribution(d2));
		//delta
		result.push_back(std::exp(-op.q*op.T)*normalDistribution(d1));
		//gamma
		result.push_back((std::exp(-op.q*op.T)) / (op.S0*op.vol*std::sqrt(op.T))*(1) / (std::sqrt(2 * pi))*std::exp(-d1*d1 / 2.0));
		//theta
		result.push_back(-(op.S0*op.vol*std::exp(-op.q*op.T)) / (2 * std::sqrt(2 * pi*op.T))*std::exp(-d1*d1 / 2.0) + op.q*op.S0*std::exp(-op.q*op.T)*normalDistribution(d1) - op.r*op.K*std::exp(-op.r*op.T)*normalDistribution(d2));
		//vega
		result.push_back(op.S0*std::exp(-op.q*op.T)*std::sqrt(op.T)*(1.0 / std::sqrt(2.0*pi))*std::exp(-d1*d1 / 2.0));

	}

	return result;
}

double BSImpliedVol(const Option& op, const double & ini1, const double & ini2, const double& marketp, const char & pc)
{
	double impvol = 0;

	std::function<double(double)> d1 = [=](double vol) {return (std::log(op.S0 / op.K) + (op.r - op.q + vol*vol / 2.0)*op.T) / (vol*std::sqrt(op.T)); };
	std::function<double(double)> d2 = [=](double vol) {return d1(vol) - vol*std::sqrt(op.T); };
	std::function<double(double)> pprice = [=](double vol) {return op.K*std::exp(-op.r*op.T)*(normalDistribution(-d2(vol))) - op.S0 *std::exp(-op.q*op.T) *(normalDistribution(-d1(vol))); };
	std::function<double(double)> cprice = [=](double vol) {return op.S0 *std::exp(-op.q*op.T) *(normalDistribution(d1(vol))) - op.K*std::exp(-op.r*op.T)*(normalDistribution(d2(vol))); };
	std::function<double(double)> psolver = [=](double vol) {return pprice(vol) - marketp; };
	std::function<double(double)> csolver = [=](double vol) {return cprice(vol) - marketp; };

	if (pc == 'p')
	{
		impvol = Secant(ini1, ini2, psolver);
	}
	else if (pc == 'c')
	{
		impvol = Secant(ini1, ini2, csolver);

	}
	else
	{
		return -1;
	}
	return impvol;
}


#endif

