#include <vector>
#include <iostream>
#include "NumberGenerators.hpp"
using namespace std;

vector<double> NumberGenerators::LinearCongruential(int num) {
	long long x = 1;
	int a = 39373;
	int k = pow(2, 31) - 1;
	int c = 0;

	vector<double> Result;
	for (int count = 0; count < num; count++) {
		x = (x*a + c) % k;
		double u = double(x) / k;
		Result.push_back(u);
	}
	return Result;
}

vector<double> NumberGenerators::InverseTransform(int num) {
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

	NumberGenerators NG;
	vector<double> UniformSamples = NG.LinearCongruential(num);

	double x = 0;
	double y = 0, r = 0;
	int cnt = 0;
	vector<double> Result;
	for (int i = 0; i < num; i++) {
		y = UniformSamples[i] - 0.5;

		if (abs(y) < 0.42) {
			r = y*y;
			x = y*(((a3*r + a2)*r + a1)*r + a0) / ((((b3*r + b2)*r + b1)*r + b0)*r + 1);
		}
		else {
			r = UniformSamples[i];
			if (y > 0) {
				r = 1 - UniformSamples[i];
			}
			r = log(-log(r));
			x = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 + r*c8)))))));
			if (y < 0) x = -x;
		}

		Result.push_back(x);
	}
	return Result;
}

vector<double> NumberGenerators::AcceptanceRejection(int num) {
	NumberGenerators Ng;
	vector<double> uniformSample = Ng.LinearCongruential(num);
	int i = 0;
	double x;
	vector<double> Results;
	while (i < num - 2) {
		x = -log(uniformSample[i]);
		if (uniformSample[i + 1] > exp(pow(x - 1, 2)*(-0.5))) {
			i = i + 1;
		}
		else {
			if (uniformSample[i + 2] <= 0.5) { x = -x; }
			Results.push_back(x);
			i = i + 1;
		}
	}
	return Results;
}

vector<double> NumberGenerators::BoxMuller(int num) {
	NumberGenerators Ng;
	vector<double> uniformSamples = Ng.LinearCongruential(num * 2);

	vector<double> Results;
	double x, y;
	double u1, u2;
	double flag = 0;
	int count = 0;

	//generate independent uniform samples
	while (count < num) {
		x = 2;
		while (x > 1) {
			u1 = 2 * uniformSamples[flag] - 1;
			u2 = 2 * uniformSamples[flag + 1] - 1;
			x = u1*u1 + u2*u2;
			flag += 2;
		}
		y = sqrt(-2 * log(x) / x);
		Results.push_back(u1*y);
		Results.push_back(u2*y);
		count += 2;
	}
	return Results;
}

double Ninverse(double x) {
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

	double result;
	double y = x - 0.5;
	if (abs(y) < 0.42) {
		double r = pow(y, 2);
		result = y*(a0 + r*(a1 + r*(a2 + r*a3))) /
			(1 + r*(b0 + r*(b1 + r*(b2 + r*b3))));
	}
	else {
		double r = x;
		if (y > 0) { r = 1 - x; };
		r = log(-log(r));
		result = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 + r*c8)))))));
		if (y < 0) { result = -result; }
	}
	return result;
}