#ifndef MONTE_CARLO
#define MONTE_CARLO

#include <iostream>
#include <Eigen\Dense>
#include "OptionPricing.hpp"
using namespace std;
using namespace Eigen;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

class MonteCarlo : public OptionPricing {
public:
	MonteCarlo() {};
	MonteCarlo(double s, double k, double t, double rr, double qq, double ss) :OptionPricing(s, k, t, rr, qq, ss) {};
	~MonteCarlo() {};

	vector<double> PathDependent_Div(int N, vector<double> Div_times, vector<double> Div_values);
};

double MeanofVector(vector<double> Input);
#endif // !MONTE_CARLO
