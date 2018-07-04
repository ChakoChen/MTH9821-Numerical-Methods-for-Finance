//BarrierOption
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>
#include "BarrierOption.hpp"
using namespace std;

void PrintVector(vector<double> Input) {
	for_each(Input.begin(), Input.end(), [](double i) {cout << i << setprecision(12) << ", "; });
	cout << endl;
}

int main() {
	double S = 42;
	double K = 40;
	double B = 35;
	double T = 7.0 / 12;
	double q = 0.03;
	double r = 0.05;
	double sigma = 0.3;
	double w = 1.2;
	double tol = pow(10, -6);

	BarrierOption Opt(S, K, B, T, q, r, sigma);
	
	double alpha = 0.4;
	int M;
	vector<double> domain;
	cout << "alpha = " << alpha << " :" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		domain = Opt.DomainDiscret(alpha, M);
		PrintVector(domain);
	}
	alpha = 4;
	cout << "alpha = " << alpha << " :" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		domain = Opt.DomainDiscret(alpha, M);
		PrintVector(domain);
	}

	vector<double> result;
	cout << "Forward Euler" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		result = Opt.ForwardEuler(0.4, M, "NO");
		PrintVector(result);

	}

	cout << "Backward Euler with alpha = 0.4" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		result = Opt.BackwardEuler_LU(0.4, M, "NO");
		PrintVector(result);
	}

	cout << "Backward Euler with alpha = 4" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		result = Opt.BackwardEuler_LU(4, M, "NO");
		PrintVector(result);
	}

	cout << "Crank Nicolson with alpha = 0.4" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		result = Opt.CrankNicolson_SOR(0.4, M, w, tol);
		PrintVector(result);
	}
	cout << "Crank Nicolson with alpha = 4" << endl;
	for (int i = 1; i <= 4; i++) {
		M = pow(4, i);
		result = Opt.CrankNicolson_SOR(4, M, w, tol);
		PrintVector(result);
	}
	
	cout << "Forward Euler with alpha = 0.4" << endl;
	M = 4;
	Opt.ForwardEuler(0.4, M, "Yes");
	cout << "Backward Euler with alpha = 0.4" << endl;
	Opt.BackwardEuler_LU(0.4, M, "Yes");


	
}