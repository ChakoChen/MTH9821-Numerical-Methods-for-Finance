/* Test_MC.cpp
 * 09/17/17
 *
 * Chapin Day
 *
 * Testing Monte Carlo simulation
 *
 * This program will generate independent uniform random numbers, then employ one of 3 different
 * methods for generating independent normal random numbers from the uniform set, and then
 * price European and American (i.e. path- and non-path-dependent options).
 *
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "Option.hpp"
#include "Black_Scholes.hpp"
#include "Monte_Carlo_Vanilla.hpp"
#include "Monte_Carlo_Path_Dependent.hpp"
#include "Option_Pricer.hpp"
#include "StopWatch.hpp"
#include <memory>
#include <cmath> // for variance
#include "Random_Number_Generator.hpp"

using pVec = std::shared_ptr<std::vector<double> >; // the return type

void print_tuple(const std::tuple<double,double,double,double,double>&);
void test_rng();
void test_vanilla_options();
void test_path_dependent_options();
void problem_3();

int main()
{
	// increase output precision
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(14);

	//test_vanilla_options();
	//test_path_dependent_options();
	problem_3();
	//test_rng();

	std::cout << "\n\nDone testing.\n";
	return 0;
}


void print_tuple(const std::tuple<double,double,double,double,double>& t)
{
	double delta,gamma,theta,vega,price;
	std::tie(delta,gamma,theta,vega,price) = t;
	std::cout << " Price: " << price << ", Delta: " << delta << ", Vega: " << vega;
}

void test_rng()
{
	////////////////// test random number generator ///////////////
	Random_Number_Generator RNG;

	unsigned int len = 10000;

	std::cout << "Generating " << len << " Uniform random numbers";
	pVec U = RNG.Uniform(len,1);
	stats(U);

	/*
	std::ofstream myfile;
	myfile.open("u.csv");
	for (auto & elem : *U)
	{
		myfile << elem << "\n";
	}
	myfile.close();
	*/

	auto iter1 = (*U).begin();
	auto iter2 = (*U).end();
	std::cout << "\nVector info: end - begin: " << iter2 - iter1 << std::endl;

	std::cout << "Generating " << len << " Normal random numbers (Inverse Transform)\n";
	pVec N1 = RNG.Inverse_Transform(len);
	stats(N1);

	std::ofstream myfile;
	myfile.open("it.csv");
	for (auto & elem : *N1)
	{
		myfile << elem << "\n";
	}
	myfile.close();

	std::cout << "Generating " << len << " Normal random numbers (Accept-Reject)\n";
	pVec N2 = RNG.Accept_Reject(len);
	stats(N2);

	std::cout << "Generating " << len << " Normal random numbers (Box-Muller):\n";
	pVec N3 = RNG.Box_Muller(len);
	stats(N3);
}

void test_vanilla_options()
{
	///////////////////////////// PLAIN VANILLA OPTIONS ///////////////////
	// option pricing parameters
	Option Call_09_42(false,1,41,42,0.75,0.01);
	Option Put_09_42(false,-1,41,42,0.75,0.01);
	std::cout << "Created options: \n" << Call_09_42 << "\n" << Put_09_42 << std::endl;
	double sigma = 0.20;
	double r = 0.03;

	StopWatch timer;

	std::shared_ptr<Option_Pricer> BS(new Black_Scholes);

	std::cout << "Call: " ;
	print_tuple(BS->Greeks(Call_09_42,sigma,r,0));
	std::cout << "Put: ";
	print_tuple(BS->Greeks(Put_09_42,sigma,r,0));

	std::shared_ptr<Option_Pricer> MC(new Monte_Carlo_Vanilla);

	int N;

	std::cout << "\n\n# Paths, Call Px, Call Delta, Call Vega, Put Px, Put Delta, Put Vega\n";

	for (int k=0; k<=9; ++k)
	{
		N = 10000 * pow(2,k);
		std::cout << N; 
		timer.StartStopWatch();
		print_tuple(MC->Greeks(Call_09_42,sigma,r,N));
		print_tuple(MC->Greeks(Put_09_42,sigma,r,N));
		timer.StopStopWatch();
		std::cout << " " << timer.GetTime() << std::endl;
	}
}


void test_path_dependent_options()
{
	// Note: option arguments: American(T/F), type (1/-1), S, K, T, q, Barrier, barrier_type
	// where barrier_type is 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in

	// down-and-out call with $35 barrier
	Option BCall(false,1,39,39,0.75,0.01,35,3);
	double sigma = 0.25;
	double r = 0.02;
	std::cout << "Created " << BCall << std::endl;

	// calculate precise value
	std::shared_ptr<Option_Pricer> BS(new Black_Scholes);
	double C_S_K = BS->Price(BCall,sigma,r);
	double a = (r - BCall.q())/(sigma*sigma) - 0.5;
	double term1 = pow(BCall.B() / BCall.S(), 2*a);
	Option tmp ( BCall);
	tmp.S(BCall.B() * BCall.B() / BCall.S()); // set price of temp option to B^2 / S
	double term2 = BS->Price(tmp,sigma,r);
	double Price_exact = C_S_K - term1 * term2;

	std::cout << "Price exact: " << Price_exact << std::endl;

	// Nk = 10000 * 2^9
	int num_paths = 50*pow(2,9); // multiplying by 200 paths below completes Nk

	// calculate monte carlo value
	std::shared_ptr<Option_Pricer> MC_Path(new Monte_Carlo_Path_Dependent);
	double Price_MC = MC_Path->Price(BCall,sigma,r,num_paths,200);
	std::cout << "MC Price: " <<  Price_MC << std::endl;
	double approx_error = std::abs(Price_exact - Price_MC);

	std::cout << "Approximation error: " << approx_error << std::endl;

	/* 
	// Just for fun, see if values of barrier calls add up to black-scholes call
	// test value of down-and-out call with down-and-in call, compared to plain call
	Option BCall2(false,1,39,39,0.75,0.01,35,4);
	std::cout << "Created " << BCall2 << std::endl;
	std::cout << "MC Price: " << MC_Path->Price(BCall2,sigma,r,num_paths,200) << std::endl;

	Option BCall3(false,1,39,39,0.75,0.01);
	std::cout << "Created" << BCall3 << std::endl;
	std::cout << "MC Price: " << MC_Path->Price(BCall3,sigma,r,num_paths,200) << std::endl;

	std::cout << "Double-checking  BS price of last call: "
		<< BS->Price(BCall3,sigma,r) << std::endl;
		*/

	std::cout << "\nProblem 3.1 - m=200; variable paths:\n";
	std::cout << "#nodes; m; n; V; |Cdao-V|; mk; nk; Vnk; |Cdao-Vnk|\n";
	int m=200;
	for (int k=0; k<10; ++k)
	{
		int n = 50*pow(2,k);
		std::cout << n*m << ";"; // print  number of nodes
		std::cout << m << ";" << n << ";"; // print m and n
		double V = MC_Path->Price(BCall,sigma,r,n,m);
		std::cout << V << ";"; // print MC option value;
		std::cout << std::abs(Price_exact - V) << ";"; // print error vs BS exact value
		double mk = ceil(pow(n*m,1.0/3)*pow(BCall.T(),2.0/3));
		double nk = floor(n*m/mk);
		std::cout << mk << ";" << nk << ";"; // print mk,nk
		V = MC_Path->Price(BCall,sigma,r,nk,mk);
		std::cout << V << ";" << std::abs(Price_exact - V) << std::endl;
	}

}

void problem_3()
{
	// Note: option arguments: American(T/F), type (1/-1), S, K, T, q, Barrier, barrier_type
	// where barrier_type is 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in

	Option Put55(false,-1,50,55,0.5,0);
	double sigma = 0.30; // volatility
	double r = 0.04; // risk-free rate

	std::shared_ptr<Option_Pricer> BS(new Black_Scholes);

	std::cout << Put55;
	std::cout << "BS price: ";
	double price_BS = BS->Price(Put55,sigma,r);
	std::cout << price_BS << std::endl;

	std::shared_ptr<Option_Pricer> MC(new Monte_Carlo_Vanilla);

	std::vector<std::string> methods;
	methods.push_back("Inverse Transform");
	methods.push_back("Accept-Reject");
	methods.push_back("Box-Muller");
	auto method = methods.begin();

	double price_MC;

	for (int i=0; i<3; ++i)
	{ // iterate through methods
		std::cout << *method++ << std::endl;
		for (int k=0; k<10; ++k)
		{ // increment number of paths
			int Nk = 10000 * pow(2,k);
			std::cout << Nk << ";";
			price_MC = MC->Price(Put55,sigma,r,Nk,i); // i chooses RNG method
			std::cout << price_MC << ";";
			std::cout << std::abs(price_BS - price_MC) << std::endl;
		}
	}
}

