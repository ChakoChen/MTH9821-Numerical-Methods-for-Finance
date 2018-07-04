//Homework 1

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "BinomialTreePricer.h"
#include "BlackScholes.h"
#include <limits>

typedef std::vector<double>(*models) (Option op, const int & N, const char & optype, const char & pc);
models functions[] = { BT,ABT,BBS,BBSR };

void question1()
{
	Option op;
	//op.K = 40;
	//op.S0 = 41;
	//op.q = 0.01;
	//op.vol = 0.3;
	//op.r = 0.03;
	//op.T = 1;

	op.K = 41;
	op.S0 = 45;
	op.q = 0.01;
	op.vol = 0.35;
	op.r = 0.025;
	op.T = 1;

	BinomialTree price(op);
	std::vector<double> bt(4);	//{price, delta, gamma, theta}
	std::vector<double> abt(4); 
	std::vector<double> bbs(4); 
	std::vector<double> bbsr(4);

	int N1 = 10;
	int N2 = 100;

	std::cout << "\nEuropean Put Option\n";
	for (int i = N1; i <= N2; i++)
	{
		bt = price.BTOptionPricer(i, 'e', 'p');
		abt = price.ABTOptionPricer(i, 'e', 'p');
		bbs = price.BBSOptionPricer(i, 'e', 'p');
		bbsr = price.BBSROptionPricer(i, 'e', 'p');
		int precision = std::numeric_limits<double>::max_digits10;
		std::cout << std::setprecision(precision) << bt[0] << "," << abt[0] << "," << bbs[0] << "," << bbsr[0] << "\n";
	}

	std::cout << "\nAmerican Put Option\n";
	for (int i = N1; i <= N2; i++)
	{
		bt = price.BTOptionPricer(i, 'a', 'p');
		abt = price.ABTOptionPricer(i, 'a', 'p');
		bbs = price.BBSOptionPricer(i, 'a', 'p');
		bbsr = price.BBSROptionPricer(i, 'a', 'p');
		int precision = std::numeric_limits<double>::max_digits10;
		std::cout << std::setprecision(precision) << bt[0] << "," << abt[0] << "," << bbs[0] << "," << bbsr[0] << "\n";
	}

	std::cout << std::endl;
}

void question2()
{
	Option op;
	op.K = 41;
	op.S0 = 45;
	op.q = 0.01;
	op.vol = 0.35;
	op.r = 0.025;
	op.T = 1;

	//Black Scholes option
	std::vector<double> bs(4);
	char pc;
	pc = 'p';
	bs = BSOption(op,pc);
	std::cout << "\nBlack Scholes Option\n";
	int precision = std::numeric_limits<double>::max_digits10;
	for_each(bs.begin(), bs.end(), [=](double x) {std::cout << std::setprecision(precision) << x << ","; });
	std::cout << "\n\n";


	//European Put Options pricer using Binomial tree, average binomial tree,
	//binomial black-scholes, and binomial black-scholes with Richardson Extrapolation
	
	//{V(N), |V(N)-V_BS|, N|V(N)-V_BS|, N^2|V(N)-V_BS|, Delta, |Delta-Delta_BS|, Gamma, |Gamma-Gamma_BS|, Theta, |Theta-Theta_BS|}
	std::vector<double> res; 
	std::vector<double> PnG(4);	//price and greeks

	//Time intervals are {10,20,40,80,160,320,640,1280}
	for (int k = 0; k < 4; k++)
	{
		for (int N = 10; N <= 1280;)
		{
			PnG = functions[k](op,N, 'e', 'p');

			//prices and their errors
			res.push_back(PnG[0]);							//V(N)
			res.push_back(std::abs(PnG[0] - bs[0]));			//|V(N)-V_BS|
			res.push_back(N*std::abs(PnG[0] - bs[0]));		//N*|V(N)-V_BS|
			res.push_back(N*N*std::abs(PnG[0] - bs[0]));		//N^2*|V(N)-V_BS|
			for (int j = 1; j < 4; j++)
			{
				res.push_back(PnG[j]);						//greeks
				res.push_back(std::abs(PnG[j] - bs[j]));		//absolute err compare to bs results
			}
			int precision = std::numeric_limits<double>::max_digits10;
			for_each(res.begin(), res.end(), [=](double x) {std::cout << std::setprecision(precision) << x << ","; });
			
			std::cout << "\n";

			res.clear();
			N = N * 2;
		}
		std::cout << "\n";
	}
}

void question3a()
{
	Option op;
	op.K = 41;
	op.S0 = 45;
	op.q = 0.01;
	op.vol = 0.35;
	op.r = 0.025;
	op.T = 1;
	BinomialTree Pricer(op);

	//Use the value of an American Put by using an average binomial tree with 10000 and 10001 steps
	std::vector<double> exact = Pricer.ABTOptionPricer(10000,'a','p');
	for_each(exact.begin(), exact.end(), [](double x) {std::cout << std::fixed << std::setprecision(6) << x << ","; });
	std::cout << "\n\n";

	//{V(N), |V(N)-V_BS|, N|V(N)-V_BS|, N^2|V(N)-V_BS|, Delta, |Delta-Delta_BS|, Gamma, |Gamma-Gamma_BS|, Theta, |Theta-Theta_BS|}
	std::vector<double> res;
	std::vector<double> PnG(4);	//price and greeks

	//Time intervals are {10,20,40,80,160,320,640,1280}
	for (int k = 0; k < 4; k++)
	{
		for (int N = 10; N <= 1280;)
		{
			PnG = functions[k](op, N, 'a', 'p');

			//prices and their errors
			res.push_back(PnG[0]);							//V(N)
			res.push_back(std::abs(PnG[0] - exact[0]));			//|V(N)-V_BS|
			res.push_back(N*std::abs(PnG[0] - exact[0]));		//N*|V(N)-V_BS|
			res.push_back(N*N*std::abs(PnG[0] - exact[0]));		//N^2*|V(N)-V_BS|
			for (int j = 1; j < 4; j++)
			{
				res.push_back(PnG[j]);						//greeks
				res.push_back(std::abs(PnG[j] - exact[j]));		//absolute err compare to bs results
			}
			int precision = std::numeric_limits<double>::max_digits10;
			for_each(res.begin(), res.end(), [=](double x) {std::cout << std::setprecision(precision) << x << ","; });

			std::cout << "\n";

			res.clear();
			N = N * 2;
		}
		std::cout << "\n";
	}
}

void question3b()
{
	Option op;
	op.K = 41;
	op.S0 = 45;
	op.q = 0.01;
	op.vol = 0.35;
	op.r = 0.025;
	op.T = 1;
	BinomialTree Pricer(op);
	char pc;
	pc = 'p';

	//Black Scholes option
	std::vector<double> bs(4);
	bs = BSOption(op,pc);

	//Use the value of an American Put by using an average binomial tree with 10000 and 10001 steps
	std::vector<double> exact = Pricer.ABTOptionPricer(10000, 'a', 'p');


	//{V(N), |V(N)-V_BS|, N|V(N)-V_BS|, N^2|V(N)-V_BS|, Delta, |Delta-Delta_BS|, Gamma, |Gamma-Gamma_BS|, Theta, |Theta-Theta_BS|}
	std::vector<double> res;
	std::vector<double> PnG(4);	//price and greeks
	std::vector<double> Euro;
	std::vector<double> Amer(4);

	//Time intervals are {10,20,40,80,160,320,640,1280}
	for (int k = 0; k < 4; k++)
	{
		for (int N = 10; N <= 1280;)
		{
			Euro = functions[k](op, N, 'e', 'p');
			Amer = functions[k](op, N, 'a', 'p');
			for (int t = 0; t < 4; t++)
			{
				PnG[t] = Amer[t] - (Euro[t] - bs[t]);
			}
			//prices and their errors
			res.push_back(PnG[0]);							//V(N)
			res.push_back(std::abs(PnG[0] - exact[0]));			//|V(N)-V_BS|
			res.push_back(N*std::abs(PnG[0] - exact[0]));		//N*|V(N)-V_BS|
			res.push_back(N*N*std::abs(PnG[0] - exact[0]));		//N^2*|V(N)-V_BS|
			for (int j = 1; j < 4; j++)
			{
				res.push_back(PnG[j]);						//greeks
				res.push_back(std::abs(PnG[j] - exact[j]));		//absolute err compare to bs results
			}
			int precision = std::numeric_limits<double>::max_digits10;
			for_each(res.begin(), res.end(), [=](double x) {std::cout << std::setprecision(precision) << x << ","; });

			std::cout << "\n";

			res.clear();
			N = N * 2;
		}
		std::cout << "\n";
	}
}

int main()
{
	std::cout << "\nHomework 1 CPP 1\n";
	question1();

	std::cout << "\n\n\nHomework 1 CPP 2\n";
	question2();

	std::cout << "\n\n\nHomework 1 CPP 3a\n";
	question3a();

	std::cout << "\n\n\nHomework 1 CPP 3b\n";
	question3b();

	system("pause");
	system("pause"); 
	system("pause");
	return 0;
}
