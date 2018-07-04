#ifndef FINITE_DIFF
#define FINITE_DIFF

#include <tuple>
#include <math.h>
#include <Eigen/Dense>
#include "OptionPricing.hpp"
#include "LinearSolver.hpp"
using namespace std;
using namespace Eigen;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

class FiniteDiff : public OptionPricing {
public:
	FiniteDiff() {};
	FiniteDiff(double s, double k, double t, double rr, double qq, double ss) :OptionPricing(s, k, t, rr, qq, ss) {};
	~FiniteDiff() {};

	tuple<vector<double>, vector<double>> DomainDiscrete_Div(double alpha, int M_1, double t_div, double alpha_2);

	tuple<vector<double>, mat, mat> ForwardEuler_Div(double alpha_1, int M_1, double t_div, double alpha_2);

	tuple<vector<double>, mat, mat> CrankNicolson_Div(double alpha_1, int M_1, double t_div, double alpha_2);
};

#endif // !FINITE_DIFF
