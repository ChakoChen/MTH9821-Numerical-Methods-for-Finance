#include <iostream>
#include "RandomNumberGenerators.hpp"

double LinearCongruentialGenerator::getUniform()
{ // get the next uniform random number

	x = (a * x + c) % k;

	return double(x) / k;	// force to return a double
};

//Linear Congruential Generator
std::vector<double> Linear_Congruential(int N) {
	std::vector<double> res(N);
	int64_t x = 1;
	int a = 39373, c = 0, k = (int)std::pow(2, 31) - 1;
	for (int i = 0; i < N; ++i) {
		x = (a * x + c) % k;
		res[i] = (double(x) / k);
	}
	return res;
}
//Acceptance-Rejection Method
std::vector<double> Acceptance_Rejection(int N) {
	std::vector<double> res;
	std::vector<double> u = Linear_Congruential(4 * N);
	for (auto it = u.begin(); it <= u.end() - 3; it = it + 3) {
		auto X = -std::log(*it);
		if (*(it + 1) > std::exp(-0.5 * (X - 1) * (X - 1)))
			continue;
		else
			if (*(it + 2) <= 0.5)
				X = -X;
		res.push_back(X);
		if (res.size() == N)
			return res;
	}
	return res;
}

//Box-Muller Method
std::vector<double> Box_Muller(int N) {

	LinearCongruentialGenerator UniformGenerator = LinearCongruentialGenerator();
	
	std::vector<double> res(N);
	
	for (auto i = 0; i < N; i = i + 2) {
	
		double U1;
		double U2;
		double X = 2;
		while (X > 1) {

			U1 = UniformGenerator.getUniform();
			U2 = UniformGenerator.getUniform();
			U1 = 2 * U1 - 1;
			U2 = 2 * U2 - 1;
			X = U1 * U1 + U2 * U2;
		
		}

		double Y = sqrt(-2 * log(X) / X);
		
		double Z1 = U1 * Y;
		
		double Z2 = U2 * Y;
		
		res[i] = Z1;
		
		res[i+1] = Z2;
	}

	return res;
}

//Inverse-Transform Method
std::vector<double> Inverse_Transform(int N) {
	std::vector<double> res;
	std::vector<double> u = Linear_Congruential(4 * N);
	double a0 = 2.50662823884,
		a1 = -18.61500062529,
		a2 = 41.39119773534,
		a3 = -25.44106049637,
		b0 = -8.47351093090,
		b1 = 23.08336743743,
		b2 = -21.06224101826,
		b3 = 3.13082909833,
		c0 = 0.3374754822726147,
		c1 = 0.9761690190917186,
		c2 = 0.1607979714918209,
		c3 = 0.0276438810333863,
		c4 = 0.0038405729373609,
		c5 = 0.0003951896511919,
		c6 = 0.0000321767881768,
		c7 = 0.0000002888167364,
		c8 = 0.0000003960315187;
	for (auto it = u.begin(); it < u.end(); ++it) {
		double y = *it - 0.5;
		double z;
		if (std::abs(y) < 0.42) {
			z = y*y;
			double numerator = ((a3 * z + a2) * z + a1) * z + a0;
			double denominator = (((b3 * z + b2) * z + b1) * z + b0) * z + 1;
			res.push_back(y * numerator / denominator);
		}
		else {
			z = *it;
			if (y > 0) {
				z = 1 - z;
			}
			z = std::log(-std::log(z));
			double temp = (((((((c7 + z * c8) * z + c6) * z + c5) * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
			if (y < 0) {
				res.push_back(-temp);
			}
			else {
				res.push_back(temp);
			}
		}
		if (res.size() == N)
			return res;
	}
	return res;
}
