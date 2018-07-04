#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include "Option.hpp"
using namespace std;

//function to print 2-dimention vector
void PrintMatrix(vector<vector<double>> Result) {
	for (int i = 0; i < Result.size(); i++) {
		for_each(Result[i].begin(), Result[i].end(), [](double i) {cout << setprecision(12)<< i << ", "; });
		cout << endl;
	}
}
//function to print 1-dimension vector
void PrintVector(vector<double> Input){
	for_each(Input.begin(), Input.end(), [](double i) {cout << setprecision(12) << i << ", "; });
}

int main() {
	double S = 42;
	double K = 40;
	double T = 0.75;
	double r = 0.04;
	double q = 0.02;
	double sigma = 0.32;
	Option Opt(S, K, T, r, q, sigma, "PUT");
	
	cout << "Backward Euler with LU and ¦Á = 0.45: " << endl;
	PriceResult BW = Opt.BackWardEuler(4, 0.45);
	PrintMatrix(BW.Umatrix);
	cout << endl;
	cout << "Crank-Nicolson with LU and ¦Á = 0.45: " << endl;
	PriceResult CN = Opt.CrankNicolson(4, 0.45);
	PrintMatrix(CN.Umatrix);
	cout << endl;

	
	cout << "\nBackward Euler with LU and ¦Á = 0.45: " << endl;
	for (int i = 1; i < 5; i++) {
		int M = pow(4, i);
		PriceResult BW = Opt.BackWardEuler(M, 0.45);
		PrintVector(BW.PointwiseConvergence);
		cout << endl;
	}
	cout << "\nCrank-Nicolson with LU and ¦Á = 0.45: " << endl;
	for (int i = 1; i < 5; i++) {
		int M = pow(4, i);
		PriceResult CN = Opt.CrankNicolson(M, 0.45);
		PrintVector(CN.PointwiseConvergence);
		cout << endl;
	}
	cout << "\nBackward Euler with LU and ¦Á = 5: " << endl;
	for (int i = 1; i < 5; i++) {
		int M = pow(4, i);
		PriceResult BW = Opt.BackWardEuler(M, 5);
		PrintVector(BW.PointwiseConvergence);
		cout << endl;
	}
	cout << "\nCrank-Nicolson with LU and ¦Á = 5: " << endl;
	for (int i = 1; i < 5; i++) {
		int M = pow(4, i);
		PriceResult CN = Opt.CrankNicolson(M, 5);
		PrintVector(CN.PointwiseConvergence);
		cout << endl;
	}

	return 0;
}