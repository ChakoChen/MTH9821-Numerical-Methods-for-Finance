#ifndef OPTION_HPP
#define OPTION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <eigen/dense>
#include "RWcsv.h"
//Option(S0,K,T,sigma,q,r)
struct Option {
	//underlying stock
	double S0;
	//Option
	double K;
	double T;
	double sigma;
	double q;
	double r;
};

class Pricer
{
private:
	Option op;
	char type;

public:
	Pricer() {}
	Pricer(Option src, char t) :op(src), type(t) {}
	~Pricer() {}

//	double PayOffPut(double St) { return std::max<double>(op.K - St, 0); }
//	double PayOffCall(double St) { return std::max<double>(St - op.K, 0); }



	double EuropreanOptionFiniteDifference(int M, double alpha_temp);

	std::vector<double> AmericanOptionFiniteDifference(int M, double alpha_temp);

	
};

double Pricer::EuropreanOptionFiniteDifference(int M, double alpha_temp)
{
	double price = 0;

	//change variables
	double a = (op.r - op.q) / (op.sigma*op.sigma) - 0.5;
	double b = (a + 1)*(a + 1) + 2 * op.q / (op.sigma*op.sigma);

	//Define the Computational Domain
	//\tao_final
	double tt_final = op.T*op.sigma*op.sigma / 2.0;

	//Domain Discretization
	double x_left = std::log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T - 3.0*op.sigma*std::sqrt(op.T);
	double x_right = std::log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T + 3.0*op.sigma*std::sqrt(op.T);
	double dt = tt_final / M;
	int N = int(std::floor((x_right - x_left) / std::sqrt(dt / alpha_temp)));
	double dx = (x_right - x_left) / N;
	double alpha = dt / (dx*dx);

	std::vector<double> x;
	std::vector<double> tt;

	//the grid
	for (int n = 0; n < N + 1; n++)
	{
		x.push_back(x_left + n*dx);
	}
	for (int m = 0; m < M + 1; m++)
	{
		tt.push_back(m*dt);
	}

	//Solve the heat equation
	Eigen::MatrixXd u(M + 1,N+1);
	
	//assign the boundary values to u
	if (type == 'p') {//put option
		for (int n = 0; n < N + 1; n++)
		{
			u(0,n) = op.K*std::exp(a*x[n])*std::max<double>(1-std::exp(x[n]),0);
		}

		for (int m = 0; m < M + 1; m++)
		{
			u(m,0) = op.K*std::exp(a*x_left+b*tt[m])*(std::exp(-2*op.r*tt[m]/(op.sigma*op.sigma))-std::exp(x_left-2*op.q*tt[m]/(op.sigma*op.sigma)));
			u(m,N) = 0;
		}

	}
	else {//call option
		for (int n = 0; n < N + 1; n++)
		{
			u(0,n) = op.K*std::exp(a*x[n])*std::max<double>(std::exp(x[n]) - 1, 0);
		}
		for (int m = 0; m < M + 1; m++)
		{
			u(0,m) = 0;
			u(m,N) = op.K*std::exp(a*x_left + b*tt[m])*(std::exp(x_left - 2 * op.q*tt[m] / (op.sigma*op.sigma))- std::exp(-2 * op.r*tt[m] / (op.sigma*op.sigma)));
		}
	}


	
	//forward finite difference
	for (int m = 0; m < M; m++)	//n
	{
		for (int n = 1; n < N; n++)	//m
		{
			u(m+1,n) = alpha*u(m,n-1) + (1 - 2 * alpha)*u(m,n) + alpha*u(m,n+1);
		}
	}

	//	Output the grid
	//save_csv<MatrixXd>(u, "european.csv");

	double x_comp = std::log(op.S0 / op.K);
	double Si = 0, Si1 = 0, Vi, Vi1;
	for (int n = 0; n < N + 1; n++)
	{
		if (x[n] <= x_comp&&x_comp < x[n + 1])
		{
			Si = op.K*std::exp(x[n]);
			Si1 = op.K*std::exp(x[n + 1]);
			Vi = std::exp(-a*x[n] - b*tt_final)*u(M,n);
			Vi1 = std::exp(-a*x[n + 1] - b*tt_final)*u(M,n + 1);
			price = ((Si1 - op.S0)*Vi + (op.S0 - Si)*Vi1) / (Si1 - Si);
			break;
		}
	}
	//output the grid

	return price;
}

std::vector<double> Pricer::AmericanOptionFiniteDifference(int M, double alpha_temp)
{
	std::vector<double> result;

	//change variables
	double a = (op.r - op.q) / (op.sigma*op.sigma) - 0.5;
	double b = (a + 1)*(a + 1) + 2 * op.q / (op.sigma*op.sigma);

	//Define the Computational Domain
	//\tao_final
	double tt_final = op.T*op.sigma*op.sigma / 2.0;

	//Domain Discretization
	double x_left = std::log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T - 3.0*op.sigma*std::sqrt(op.T);
	double x_right = std::log(op.S0 / op.K) + (op.r - op.q - op.sigma*op.sigma / 2.0)*op.T + 3.0*op.sigma*std::sqrt(op.T);
	double dt = tt_final / M;
	int N = int(std::floor((x_right - x_left) / std::sqrt(dt / alpha_temp)));
	double dx = (x_right - x_left) / N;
	double alpha = dt / (dx*dx);

	std::vector<double> x;
	std::vector<double> tt;

	//the grid
	for (int n = 0; n < N + 1; n++)
	{
		x.push_back(x_left + n*dx);
	}
	for (int m = 0; m < M + 1; m++)
	{
		tt.push_back(m*dt);
	}

	//Solve the heat equation
	Eigen::MatrixXd u(M + 1, N + 1);

	//assign the boundary values to u
	if (type == 'p') {//put option
		for (int n = 0; n < N + 1; n++)
		{
			u(0, n) = op.K*std::exp(a*x[n])*std::max<double>(1 - std::exp(x[n]), 0);
		}

		for (int m = 0; m < M + 1; m++)
		{
			u(m, 0) = op.K*std::exp(a*x_left + b*tt[m])*(1-std::exp(x_left));
			u(m, N) = 0;
		}

	}
	else {//call option
		for (int n = 0; n < N + 1; n++)
		{
			u(0, n) = op.K*std::exp(a*x[n])*std::max<double>(std::exp(x[n]) - 1, 0);
		}
		for (int m = 0; m < M + 1; m++)
		{
			u(0, m) = 0;
			u(m, N) = op.K*std::exp(a*x_left + b*tt[m])*(std::exp(x_left) - 1);
		}
	}

	Eigen::MatrixXd early_ex_premium(M+1,N+1);

	for (int m = 1; m < M+1; m++)
	{
		for (int n = 1; n < N; n++)
		{
			early_ex_premium(m, n) = op.K*std::exp(a*x[n] + b*tt[m])*std::max<double>(1 - std::exp(x[n]), 0);
		}
	}

	//forward finite difference
	for (int m = 0; m < M; m++)	//n
	{
		for (int n = 1; n < N; n++)	//m
		{
			u(m + 1, n) = std::max<double>(alpha*u(m, n - 1) + (1 - 2 * alpha)*u(m, n) + alpha*u(m, n + 1), early_ex_premium(m+1, n));
		}
	}

	//	Output the grid
	save_csv<MatrixXd>(u, "M4American.csv");

	double x_comp = std::log(op.S0 / op.K);

	double price1 = 0, price2 = 0;
	int ps = 0;

	for (int n = 0; n < N + 1; n++)
	{
		if (x[n] <= x_comp&&x_comp < x[n + 1])
		{
			ps = n;
			break;
		}
	}
	
	//linear interpolate to approximate the price
	double Si = op.K*std::exp(x[ps]);
	double Si1 = op.K*std::exp(x[ps + 1]);
	double Vi = std::exp(-a*x[ps] - b*tt_final)*u(M, ps);
	double Vi1 = std::exp(-a*x[ps + 1] - b*tt_final)*u(M, ps + 1);
	price1 = ((Si1 - op.S0)*Vi + (op.S0 - Si)*Vi1) / (Si1 - Si);
	result.push_back(price1);

	//another way
	double ui = ((x[ps + 1] - x_comp)*u(M, ps) + (x_comp - x[ps])*u(M, ps + 1)) / (x[ps + 1] - x[ps]);
	price2 = std::exp(-a*x_comp - b*tt_final)*ui;
	result.push_back(price2);

	//Greeks
	double Siminus1 = op.K*std::exp(x[ps - 1]);
	double Si2 = op.K*std::exp(x[ps + 2]);
	double Viminus1 = std::exp(-a*x[ps - 1] - b*tt_final)*u(M, ps - 1);
	double Vi2 = std::exp(-a*x[ps + 2] - b*tt_final)*u(M, ps + 2);

	double delta = (Vi1 - Vi) / (Si1 - Si);
	double gamma = ((Vi2 - Vi1) / (Si2 - Si1) - (Vi - Viminus1) / (Si - Siminus1)) / ((Si2 + Si1) / 2.0 - (Si + Siminus1) / 2.0);

	double Vidt = std::exp(-a*x[ps] - b*(tt_final - dt))*u(M - 1, ps);
	double Vi1dt = std::exp(-a*x[ps + 1] - b*(tt_final - dt))*u(M - 1, ps + 1);
	double Vappdt = ((Si1 - op.S0)*Vidt + (op.S0 - Si)*Vi1dt) / (Si1 - Si);
	double theta = (Vappdt - price1) / (2.0*dt / (op.sigma*op.sigma));

	result.push_back(delta);
	result.push_back(gamma);
	result.push_back(theta);

	if (M == 16)
	{
		Eigen::MatrixXd opt(17, 1);
		int Nopt = 0;
		double Sopt = 0;
		for(int m=0;m<M+1;m++)
		{
			for (int n = 0; n < N + 1; n++)
			{
				if (std::abs(u(m, n) - early_ex_premium(m, n)) < 0.000001 && u(m,n + 1) > early_ex_premium(m,n + 1))
				{
					Nopt = n;
					break;
				}
			}
			Sopt = (op.K*std::exp(x[Nopt]) + op.K*std::exp(x[Nopt + 1])) / 2.0;
			opt(m, 0) = Sopt;
		}
		save_csv(opt, "OptimalDoman.csv");
	}

	return result;

}
#endif
