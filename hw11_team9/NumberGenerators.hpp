#ifndef RANDOMGENERATORS_HPP
#define RANDOMGENERATORS_HPP
#include <vector>
using namespace std;

class NumberGenerators {
public:
	NumberGenerators() {};
	~NumberGenerators() {};

	vector<double> LinearCongruential(int num);
	vector<double> InverseTransform(int num);
	vector<double> AcceptanceRejection(int num);
	vector<double> BoxMuller(int num);
};
double Ninverse(double x);


#endif // !RANDOMGENERATORS_HPP