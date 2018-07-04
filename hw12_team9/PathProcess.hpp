#include <tuple>
#include <algorithm>
#include <memory>
#include <string>
#include <functional>
#include <vector>
#include <iostream>
#include "AmericanPut.hpp"
#include "LinearSolver.hpp"
#include "RandomNumberGenerators.hpp"
#include <boost/math/special_functions/laguerre.hpp>

double PathProcess(std::shared_ptr<AmericanPut> Optr, int M, int N, int B, std::string random_generator, std::string method)
{
	
	double r = Optr->Interest();

	double q = Optr->Dividend();

	double sig = Optr->Volatility();

	double T = Optr->Maturity();

	double K = Optr->Strike();

	double S0 = Optr->SpotPrice();

	double delta_t = T / M;

	std::vector<double> z;

	vec C = vec::Zero(N);

	vec F = vec::Zero(N);

	vec E = vec::Zero(N);


	if (random_generator == "Inverse Transform")

		z = Inverse_Transform(N * M);

	else if (random_generator == "Acceptance Rejection")

		z = Acceptance_Rejection(N * M);

	else if (random_generator == "Box Muller") {

		z = Box_Muller(N * M);

	}

	std::vector<std::vector<double>> S(M + 1, std::vector<double>(N));

	for (int j = 0; j < N; ++j) {

		S[0][j] = S0;
	
	}

	std::vector<double> runningsum(N);

	for (int i = 1; i <= M; ++i) {

		for (int j = 0; j < N; ++j) {

			runningsum[j] += z[j * M + (i-1)];

			S[i][j] = S0 * std::exp((r - sig * sig / 2) * i * delta_t + sig * std::sqrt(delta_t) * runningsum[j]);
		
		}

	}

	for (int j = 0; j < N; ++j) {
		
		F(j) = C(j) = std::max<double>(0.0, K - S[M][j]);
	
	}

	mat basis = mat::Zero(N, B);

	for (int i = M - 1; i >= 1; --i) {

		for (int j = 0; j < N; ++j) {

			for (int k = 0; k < B; ++k) {

				if (method == "Power Basis")

					basis(j, k) = std::pow(S[i][j], k);
				
				if (method == "Laguerre Polynomials")

					basis(j, k) = std::exp(- S[i][j]/2) * boost::math::laguerre(k, S[i][j]);

			}

		}

		for (int j = 0; j < N; ++j) {

			C(j) = std::max<double>(0.0, K - S[i][j]);

		}

		vec sol = basis.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(F * std::exp(-delta_t * r));

		//vec sol = regressor_solver(basis, F * std::exp(-delta_t * r));

		E = basis * sol;

		for (int j = 0; j < N; ++j) {

			if (C(j) >= E(j))

				F(j) = C(j);

			else

				F(j) = F(j) * std::exp(-delta_t * r);
		}

	}

	return std::max<double>(E(0), C(0));
}
