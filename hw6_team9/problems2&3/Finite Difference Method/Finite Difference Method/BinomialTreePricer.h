#ifndef BINOMIALTREEPRICER_HPP
#define BINOMIALTREEPRICER_HPP



#include <iostream>
#include <cmath>
#include <vector>
#include "BlackScholes.h"

class BinomialTree
{
private:
	Option op;
public:
	BinomialTree() {}
	BinomialTree(double strike,double spot,double dividend,double volatility,double rate,double maturity) {
		op.K = strike;
		op.S0 = spot;
		op.q = dividend;
		op.sigma = volatility;
		op.r = rate;
		op.T = maturity;
	}

	BinomialTree(Option src):op(src) {}

	~BinomialTree() {}

	double SpotPrice(double u, double d, int up, int down);

	std::vector<double> BTOptionPricer(const int & N, const char & optype, const char & pc);
	std::vector<double> ABTOptionPricer(const int & N, const char & optype, const char & pc);
	std::vector<double> BBSOptionPricer(const int & N, const char & optype, const char & pc);
	std::vector<double> BBSROptionPricer(const int & N, const char & optype, const char & pc);

	
};

double BinomialTree::SpotPrice(double u,double d,int up, int down)
{
	return op.S0*std::pow(u, up)*std::pow(d, down);
}

//Binomial Tree Option pricer
//optype defines whether the option is European (e) or American (a); pc defines whether the option is call (c) or put (p)
std::vector<double> BinomialTree::BTOptionPricer(const int & N, const char & optype, const char & pc)
{
	double deltaT = op.T / (double)N;
	double u = std::exp(op.sigma*std::sqrt(deltaT));
	double d = 1.0 / u;
	double p1 = (std::exp((op.r-op.q)*deltaT) - d) / (u - d);
	double p2 = (u - std::exp((op.r-op.q)*deltaT)) / (u - d);

	std::vector<double> val1;
	std::vector<double> val2;
	double tmp1 = 0, tmp2 = 0, tmp3 = 0;
	double delta = 0, gamma = 0, theta = 0,Vud=0;
	//Payoff at time N
	for (int i = 0; i < N+1; ++i)
	{
		tmp1 = SpotPrice(u, d, N - i, i);
		if(pc=='p'){
			tmp2 = (tmp1 < op.K) ? (op.K - tmp1) : 0;	//put payoff
		}
		else {
			tmp2 = (tmp1 > op.K) ? (tmp1- op.K) : 0;	//call payoff
		}
		val1.push_back(tmp2);
	}
	for (int k = N - 1; k >= 0; --k)
	{
		for (int i = 0; i < k+1; ++i)
		{
			

			val2.resize(k+1);
			if(optype=='e'){//European option
				tmp2 = std::exp(-op.r*deltaT)*(p1*val1[i] + p2*val1[i + 1]);
				val2[i] = tmp2;
			}
			else{//American option
				tmp1 = SpotPrice(u, d, k - i, i);
				tmp2 = std::exp(-op.r*deltaT)*(p1*val1[i] + p2*val1[i + 1]);
				if (pc == 'p') {
					tmp3 = (op.K - tmp1 > tmp2) ? (op.K - tmp1) : tmp2;
				}
				else {
					tmp3 = (tmp1- op.K> tmp2) ? (tmp1 - op.K) : tmp2;
				}
				val2[i] = tmp3;
			}
		}
		if (k == 1)
		{
			delta = (val2[0] - val2[1]) / (op.S0*u - op.S0*d);
		}
		if (k == 2)
		{
			gamma = ((val2[0] - val2[1]) / (op.S0*u*u - op.S0*u*d) - (val2[1] - val2[2]) / (op.S0*u*d - op.S0*d*d)) / ((op.S0*u*u - op.S0*d*d) / 2.0);
			Vud = val2[1];
		}
		val1.resize(k + 1);
		val1 = val2;
	}
	theta = (Vud - val1[0]) / (2.0*deltaT);

	std::vector<double> ret;
	ret.push_back(val1[0]);	//price
	ret.push_back(delta);	//delta
	ret.push_back(gamma);	//gamma
	ret.push_back(theta);	//theta
	return ret;
}

//Average binomial tree model option pricer
std::vector<double> BinomialTree::ABTOptionPricer(const int & N, const char & optype, const char & pc)
{
	std::vector<double> oppN(4);
	std::vector<double> oppNI(4);
	oppN = BTOptionPricer(N, optype, pc);
	oppNI = BTOptionPricer(N+1, optype, pc);
	std::vector<double> result(4);
	for (int i= 0; i < 4; i++)
		result[i] = (oppN[i] + oppNI[i]) / 2.0;
	return result;
}

//European Option Price calculated using Binomial Black Scholes model
std::vector<double> BinomialTree::BBSOptionPricer(const int & N, const char & optype, const char & pc)
{

	double deltaT = op.T / (double)N;
	double u = std::exp(op.sigma*std::sqrt(deltaT));
	double d = 1.0 / u;
	double p1 = (std::exp((op.r - op.q)*deltaT) - d) / (u - d);
	double p2 = (u - std::exp((op.r - op.q)*deltaT)) / (u - d);

	std::vector<double> P(4);

	std::vector<double> val1;
	std::vector<double> val2;
	double tmp1 = 0, tmp2 = 0, tmp3 = 0;

	Option opBS=op;	//created for black scholes pricer

	double delta = 0, gamma = 0, theta = 0, Vud = 0;

	//Payoff at time N-1
	for (int i = 0; i < N; ++i)
	{
		opBS.S0 = SpotPrice(u, d, N - 1 - i, i);
		opBS.T = deltaT;
		P = BSOption(opBS, pc);
		if (optype == 'a')
		{//Payoff at time N-1
			if (pc == 'p')
			{
				tmp1 = (opBS.K - opBS.S0 > 0) ? opBS.K - opBS.S0 : 0;
			}
			else
			{
				tmp1 = (opBS.S0- opBS.K > 0) ? opBS.S0 - opBS.K : 0;
			}
			val1.push_back(P[0] > tmp1 ? P[0] : tmp1);

		}
		else
		{
			val1.push_back(P[0]);
		}
	}
	//Option Pricer
	for (int k = N - 2; k >= 0; --k)
	{
		for (int i = 0; i < k + 1; ++i)
		{
			val2.resize(k + 1);
			if (optype == 'e') {//European Option
				tmp3 = std::exp(-op.r*deltaT)*(p1*val1[i] + p2*val1[i + 1]);
				val2[i] = tmp3;
			}
			else {//American Option
				tmp1 = SpotPrice(u, d, k - i, i);
				tmp2 = std::exp(-op.r*deltaT)*(p1*val1[i] + p2*val1[i + 1]);
				if(pc=='p'){
					tmp3 = (op.K - tmp1 > 0) ? (op.K - tmp1) : 0;
				}
				else {
					tmp3 = (tmp1- op.K> 0) ? (tmp1 - op.K) : 0;
				}
				val2[i] = (tmp3 > tmp2) ? tmp3 : tmp2;
			}

		}
		if (k == 1)
		{
			delta = (val2[0] - val2[1]) / (op.S0*u - op.S0*d);
		}
		if (k == 2)
		{
			gamma = ((val2[0] - val2[1]) / (op.S0*u*u - op.S0*u*d) - (val2[1] - val2[2]) / (op.S0*u*d - op.S0*d*d)) / ((op.S0*u*u - op.S0*d*d) / 2.0);
			Vud = val2[1];
		}

		val1.resize(k + 1);
		val1 = val2;
	}
	theta = (Vud - val1[0]) / (2.0*deltaT);

	std::vector<double> ret;
	ret.push_back(val1[0]);
	ret.push_back(delta);
	ret.push_back(gamma);
	ret.push_back(theta);
	return ret;
}

//BBSR model
std::vector<double> BinomialTree::BBSROptionPricer(const int & N, const char & optype, const char & pc)
{
	std::vector<double> oppN(4);
	std::vector<double> oppNH(4);
	oppN = BBSOptionPricer(N, optype, pc);
	int half = 0;
	if (N % 2 == 0) { 
		half = N / 2; 
	}
	else {
		half = (N + 1) / 2;
	}
	oppNH = BBSOptionPricer(half, optype, pc);
	std::vector<double> result(4);
	for (int i = 0; i < 4; i++)
	{
		result[i] = 2.0*oppN[i] - oppNH[i];
	}
	return result;
}

//Create seperate functions in order to have an array of function pointers
std::vector<double> BT( Option op, const int & N, const char & optype, const char & pc)
{
	BinomialTree pricer(op);
	return pricer.BTOptionPricer(N, optype, pc);
}
std::vector<double> ABT( Option op, const int & N, const char & optype, const char & pc)
{
	BinomialTree pricer(op);
	return pricer.ABTOptionPricer(N, optype, pc);
}
std::vector<double> BBS( Option op, const int & N, const char & optype, const char & pc)
{
	BinomialTree pricer(op);
	return pricer.BBSOptionPricer(N, optype, pc);
}
std::vector<double> BBSR(Option op, const int & N, const char & optype, const char & pc)
{
	BinomialTree pricer(op);
	return pricer.BBSROptionPricer(N, optype, pc);
}


#endif
