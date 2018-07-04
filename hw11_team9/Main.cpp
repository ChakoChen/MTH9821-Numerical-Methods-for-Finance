#include <iostream>
#include <iomanip>
#include <math.h>
#include <tuple>
#include <string>
#include <Eigen\Dense>
#include "OptionPricing.hpp"
#include "BinomialPricing.hpp"
#include "FiniteDiff.hpp"
#include "MonteCarloPricing.hpp"
#include "NumberGenerators.hpp"
#include "LinearSolver.hpp"
using namespace std;
using namespace Eigen;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;


void PrintVector(vector<double> Input) {
	for_each(Input.begin(), Input.end(), [](double i) {cout << setprecision(14) << i << " "; });
	cout << endl;
}

void BTTest() {
	double S = 50;
	double K = 55.55;
	double sigma = 0.3;
	double r = 0.02;
	double q = 0.01;
	vector<double> times = { 2.0 / 12, 4.0 / 12, 6.0 / 12 };
	vector<double> values = { 0.5, 0.01, 0.75 };

	vector<double> OLD, NEW;
	int N;
	double tol;


	
	//Part 1
	cout << "Part 1: " << endl;
	Binomial SOpt(S, K, 3.0 / 12, r, q, sigma);
	Binomial LOpt(S, K, 7.0 / 12, r, q, sigma);
	
	//European Prop 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.European_PropDiv("PUT", N);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.European_PropDiv("PUT", N);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "European: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	//American Prop 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.American_PropDiv("PUT", N, 2.0/12);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.American_PropDiv("PUT", N, 2.0/12);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	
	
	//European Prop 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.European_PropDiv_Mult("PUT", N, 3);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW =LOpt.European_PropDiv_Mult("PUT", N, 3);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "European: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	//American  Prop 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.American_PropDiv_Mult("PUT", N, times);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.American_PropDiv_Mult("PUT", N, times);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	
	
	
	//Part2
	cout << "Part 2: " << endl;
	double v_div = 0.5;
	
	//European Fix 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.European_FixDiv("PUT", N, 2., v_div);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.European_FixDiv("PUT", N, 2.0/12, v_div);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "European: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	//American Fix 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.American_FixDiv("PUT", N, 2.0/12, v_div);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.American_FixDiv("PUT", N, 2.0/12, v_div);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "American: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	
	//European Fix 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.European_FixDiv_Mult("PUT", N, times, 0.5);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.European_FixDiv_Mult("PUT", N, times, 0.5);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "European: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	//American Fix 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.American_FixDiv_Mult("PUT", N, times, 0.5);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.American_FixDiv_Mult("PUT", N, times, 0.5);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;

	
	//part 3

	//European Complex 7M
	cout << "Part 3:" << endl;
	N = 7;
	tol = 10.1;
	OLD = LOpt.European_Complex("PUT", N, times, values);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.European_Complex("PUT", N, times, values);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "European: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	//American Complex 7M
	N = 7;
	tol = 10.1;
	
	OLD = LOpt.American_Complex("PUT", N, times, values);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.American_Complex("PUT", N, times, values);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(14) << OLD[0] << " " << N << " " << setprecision(14) << OLD[1]
		<< " " << setprecision(14) << OLD[2] << " " << setprecision(14) << OLD[3] << endl;
	
}

void FDTest() {
	//Finite Difference
	double S = 52;
	double sigma = 0.3;
	double K = 50;
	double T = 1;
	double r = 0.03;
	double q = 0.02;
	double t_div = 5.0 / 12;

	FiniteDiff Opt(S, K, T, r, q, sigma);

	double alpha = 0.4;
	int M_1 = 4;
	auto res = Opt.ForwardEuler_Div(alpha, M_1, t_div, alpha);

	cout << "In the interval [0, tao_div): " << endl;
	cout << setprecision(14) << get<1>(res) << endl;
	cout << "In the interval [tao_div, tao_final]:" << endl;
	cout << setprecision(14) << get<2>(res) << endl;

	cout << "\nForward Euler with alpha = 0.4: " << endl;
	for (int i = 1; i <= 4; i++) {
		M_1 = pow(4, i);
		auto res = Opt.ForwardEuler_Div(alpha, M_1, t_div, alpha);
		PrintVector(get<0>(res));
	}

	cout << "\nCrankNicolson with alpha = 4: " << endl;
	alpha = 4.0;
	for (int i = 1; i <= 4; i++) {
		M_1 = pow(4, i);
		auto res = Opt.CrankNicolson_Div(alpha, M_1, t_div, alpha);
		PrintVector(get<0>(res));
	}

	cout << "Domain with alpha = 0.4: " << endl;
	alpha = 0.4;
	vector<double> a, b, d;
	for (int i = 1; i <= 4; i++) {
		M_1 = pow(4, i);
		auto domain = Opt.DomainDiscrete_Div(alpha, M_1, t_div, alpha);
		a = get<0>(domain);
		b = get<1>(domain);
		d = { b[7], b[0],a[3] + a[4], a[1],a[2], b[1],b[2], 0.5*sigma*sigma*(T - t_div),a[6], b[6], a[5] };
		PrintVector(d);
	}
	cout << "Domain with alpha = 4: " << endl;
	alpha = 4;
	for (int i = 1; i <= 4; i++) {
		M_1 = pow(4, i);
		auto domain = Opt.DomainDiscrete_Div(alpha, M_1, t_div, alpha);
		a = get<0>(domain);
		b = get<1>(domain);
		d = { b[7], b[0],a[3] + a[4], a[1],a[2], b[1],b[2], 0.5*sigma*sigma*(T - t_div),a[6], b[6], a[5] };
		PrintVector(d);
	}


}

void MCTest() {
	//Monte Carlo
	double S = 50;
	double K = 55.55;
	double T = 7.0 / 12;
	double sigma = 0.3;
	double r = 0.02;

	vector<double> div_times = { 2.0 / 12, 4.0 / 12, 6.0 / 12 };
	vector<double> div_values = { 0.5, 0.01, 0.75 };

	int N = 4 * 10000;

	MonteCarlo Opt(S, K, T, r, 0.0, sigma);
	vector<double> res;
	
	for (int k = 0; k <= 9; k++) {
		N = 4 * 10000 * pow(2, k);
		res = Opt.PathDependent_Div(N, div_times, div_values);
		PrintVector(res);
	}

}

int main() {
	cout << "Binomial Tree: " << endl;
	BTTest();
	
	cout << "Monte Carlo: " << endl;
	MCTest();

	cout << "Finite Difference: " << endl;
	FDTest();
	
	return 0;

}