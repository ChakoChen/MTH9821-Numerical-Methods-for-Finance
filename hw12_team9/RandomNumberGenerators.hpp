#include <string>
#include <vector>
#include <cmath>

class LinearCongruentialGenerator
{
	// class for Linear Congruential Generator
	// x_i+1 = a * x_i + c (mod k)
private:
	const int a = 39373;
	const int c = 0;
	const unsigned int k = (unsigned int)std::pow(2, 31) - 1;
	int64_t x = 1;

public:
	LinearCongruentialGenerator() {};
	virtual ~LinearCongruentialGenerator() {};
	double getUniform();

};

//Linear Congruential Generator
std::vector<double> Linear_Congruential(int);

//Acceptance-Rejection Method
std::vector<double> Acceptance_Rejection(int);

//Box-Muller Method
std::vector<double> Box_Muller(int);

//Inverse-Transform Method
std::vector<double> Inverse_Transform(int);

