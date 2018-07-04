#ifndef RandomNumberGenerator_HPP
#define RandomNumberGenerator_HPP

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>

//class NormalGenerator
//{
//private:
//	std::vector<double> unif;
//	std::vector<double> normalit;
//	std::vector<double> normalar;
//	std::vector<double> normalbm;
//
//public:
//	NormalGenerator() { 
//		unif.clear(); 
//		normalit.clear(); 
//		normalar.clear();
//		normalbm.clear();
//	}
//	~NormalGenerator() {}
//
//	std::vector<double> getunif() { return unif; }
//	std::vector<double> getnormal(std::string s) { 
//		if (s == "it")
//			return normalit;
//		else if (s == "ar")
//			return normalar;
//		else if (s == "bm")
//			return normalbm;
//	}
//
//	void LinearCongruential(int s, int n);
//	void InverseTransform(int n);
//	void AcceptanceRejection(int n);
//	void BoxMuller(int n);
//
//};

std::tuple<double,int> LinearCongruential(int s)
{//generate independent samples from uniform distribution on [0,1]
	long long x = s;	//seed
	double u = 0;	//random number
	double seed = 0;

	int k = pow(2, 31) - 1;
	int a = 39373;
	int c = 0;

	x = (a*x + c) % k;
	u = (double)x / k;

	std::tuple<double,int> unif = { u,x };
	return unif;
}

std::vector<double> NLinearCongruential(int s,int n)
{
	std::vector<double> unif;
	std::tuple<double, int> UnifSeed;

	int cnt = 0;

	for (int i = 0; i < n; i++)
	{
		UnifSeed = LinearCongruential(s);
		s = std::get<1>(UnifSeed);
		unif.push_back(std::get<0>(UnifSeed));

		cnt++;

		//for degugging, print # of uniform samples that have been generated
		if (i == 10000 * pow(2, cnt)-1)
		{
			std::cout << i+1 << " random uniforms generated" << std::endl;
			cnt++;
		}
	}
	return unif;
}

std::vector<double> InverseTransform(int s,int n)
{//Use the result of uniform distribution samples to generate normal distribution samples
	//constants
	double a0 = 2.50662823884;
	double a1 = -18.61500062529;
	double a2 = 41.39119773534;
	double a3 = -25.44106049637;
	double b0 = -8.47351093090;
	double b1 = 23.08336743743;
	double b2 = -21.06224101826;
	double b3 = 3.13082909833;
	double c0 = 0.337475482276147;
	double c1 = 0.9761690190917186;
	double c2 = 0.1607979714918209;
	double c3 = 0.0276438810333863;
	double c4 = 0.0038405729373609;
	double c5 = 0.0003951896511919;
	double c6 = 0.0000321767881768;
	double c7 = 0.000002888167364;
	double c8 = 0.0000003960315187;


	double x = 0;
	double y = 0, r = 0;
	int cnt = 0;

	//generate uniform random variables
	std::vector<double> unif = NLinearCongruential(s, n);
	std::vector<double> normalit;
		
	for(int i=0;i<n;i++)
	{
		y = unif[i] - 0.5;

		if (abs(y) < 0.42)
		{
			r = y*y;
			x = y*(((a3*r + a2)*r + a1)*r + a0) / ((((b3*r + b2)*r + b1)*r + b0)*r + 1);
		}
		else
		{
			r = unif[i];
			if (y > 0)
			{
				r = 1 - unif[i];
			}
			r = log(-log(r));
			x = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 + r*c8)))))));
			if (y < 0) x = -x;
		}

		normalit.push_back(x);

		//for debugging
		if (i == 10000 * pow(2, cnt)-1)
		{
			std::cout << i+1 << " random normals generated" << std::endl;
			cnt++;
		}

	}
	return normalit;
}

std::vector<double> AcceptanceRejection(int s,int n)
{
	//generate uniform distribution samples
	double x = 0;
	int i = 0;
	std::tuple<double,int> UnifSeed;
	std::vector<double> unif(3);

	std::vector<double> normalar;

	for(int cnt=0;cnt<n;)
	{
		for (int j = 0; j < 3; j++)
		{
			UnifSeed = LinearCongruential(s);
			s = std::get<1>(UnifSeed);
			unif[j]=std::get<0>(UnifSeed);
		}
		
		x = -log(unif[0]);

		if (unif[1] > exp(-0.5*(x - 1.0)*(x - 1.0)))
		{
			continue;
		}
		else if (unif[2] <= 0.5)
		{
			x = -x;
		}
		normalar.push_back(x);
		cnt++;

	}
	return normalar;
}

std::vector<double> BoxMuller(int s,int n)
{
	//generate independent uniform samples
	std::tuple<double, int> UnifSeed;
	std::vector<double> unif(2);

	double u1 = 0, u2 = 0;
	double x = 0, y = 0;
	double z1 = 0, z2 = 0;

	std::vector<double> normalbm;

	for (int cnt = 0; cnt < n;)
	{
		for (int j = 0; j < 2; j++)
		{
			UnifSeed = LinearCongruential(s);
			s = std::get<1>(UnifSeed);
			unif[j]=std::get<0>(UnifSeed);
		}
		u1 = unif[0];
		u2 = unif[1];
		u1 = 2.0 * u1 - 1.0;
		u2 = 2.0 * u2 - 1.0;
		x = u1*u1 + u2*u2;
		if (x <= 1)
		{
			y = sqrt(-2.0*log(x) / x);
			z1 = u1*y;
			z2 = u2*y;
			normalbm.push_back(z1);
			normalbm.push_back(z2);
			cnt = cnt + 2;


		}
	}
	return normalbm;
}
#endif
