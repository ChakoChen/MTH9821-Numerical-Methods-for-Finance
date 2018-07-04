/* Monte_Carlo_Path_Dependent.cpp
 * 09/17/17
 *
 * Chapin Day
 *
 * This will return the value of an option that depends on a lognormally-distributed
 * underlying asset.
 *
 *
 */

#include "Monte_Carlo_Path_Dependent.hpp"
#include <fstream>

using otuple = std::tuple<double,double,double,double,double> ;

double Monte_Carlo_Path_Dependent::Price (Option& o, double sigma, double r, int num_steps, int param2) 
{
	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();
	double B = o.B(); // option Barrier
	int N = num_steps; // number of paths
	int m = param2; // number of intervals
	double dt = T / m;
	double rad_dt = sqrt(dt);

	std::function<bool(double,double)> barrier_test; // function to compute barrier result for path
	// argument: min,max for path; returns 1 or 0

	// define Option Barrier tests
	// 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in
	switch(o.barrier_type())
	{
		case 1: // up-and-out
			barrier_test = [&B] (double min,double max)->int { return (max>B) ? false : true; };
			break;
		case 2: // up-and-in
			barrier_test = [&B] (double min, double max) { return (max>B) ? true : false; };
			break;
		case 3: // down-and-out
			barrier_test = [&B] (double min, double max) { return (min<B) ? false : true; };
			break;
		case 4: // down-and-in
			barrier_test = [&B] (double min, double max) { return (min<B) ? true : false; };
			break;
		default: // no barrier test, include path
			barrier_test = [] (double min, double max) { return true; } ;
	}


	Random_Number_Generator RNG;
	auto Z = RNG.Inverse_Transform(N*m); // generate normal vector
	auto z_iter = (*Z).begin(); // points to current normal random variable

	std::shared_ptr<std::vector<double> > result(new std::vector<double>); // for option values

	// initialize option calculation variables
	double price{0}; // option value
	double S_prev = S; // previous security value begins with starting value	
	double S_curr; // current security value
	double S_max{0},S_min{99999}; // maximum, minimum security values along path

	//std::cout << "S_prev: " << S_prev << ", dt: " << dt << ", rad_dt: " << rad_dt << ", m: " << m << std::endl;
	std::shared_ptr<std::vector<std::tuple<int,double,double,double,double> > > check_vector(new std::vector<std::tuple<int,double,double,double,double> >);

	for (int path=0; path<N; ++path)
	{
		for (int interval=0; interval<m; ++interval)
		{ // each path
			// calculate current security price, node by node
			S_curr = S_prev*pow(E,(r-q-sigma*sigma/2)*dt + sigma*rad_dt * *z_iter++);
			// remember minimum and maximum along path for barrier test
			S_min = (S_curr < S_min) ? S_curr : S_min;
			S_max = (S_curr > S_max) ? S_curr : S_max;
			S_prev = S_curr;
		}
		if (barrier_test(S_min,S_max)) // test whether path went through barrier
		{ // option value to be included
			result->push_back(pow(E,-r*T)*o.Payoff(S_curr));
			check_vector->push_back(std::make_tuple(path,S_curr,pow(E,-r*T)*o.Payoff(S_curr),S_min,S_max));
		}
		else
		{ // option failed barrier test; exclude value
			result->push_back(0); 
			//check_vector->push_back(std::make_tuple(path,S_curr,0.0,S_min,S_max));
		}
		//std::cout << "Si: " << S_curr << "; max: " << S_max << ", min: " << S_min
		//	<< "; test: " << barrier_test(S_min,S_max) << std::endl;
		S_prev = S; // start new path from original security price
		S_max = 0; // reset maximum
		S_min = 99999; // reset minimum
	}

	/*
	// write check_vector
	std::ofstream myfile;
	myfile.open("foo.csv");
	myfile << "path;S_curr;o.Payoff;S_min;S_max\n";
	for (auto & elem : *check_vector)
	{
		int p1;
		double a1,a2,a3,a4;
		std::tie(p1,a1,a2,a3,a4) = elem;
		myfile << p1 << ";" << a1 << ";" << a2 << ";" << a3 << ";" << a4 << "\n";
	}
	myfile.close();
	*/


	// calculate average price
	for (auto & elem : *result)
	{
		price += elem;
	}

	/*
	// write another check_vector
	myfile.open("foo1.csv");
	myfile << "count;price\n";
	int counter=0;
	for (auto & elem : *result)
	{
		myfile << counter++ << ";" << elem << "\n";
	}
	myfile.close();
	*/

	return price/N;
}


otuple Monte_Carlo_Path_Dependent::Greeks (Option& o, double sigma, double r, int num_steps, int p2) 
{ // delta,gamma,theta,vega

   return std::make_tuple(0,0,0,0,0); // Test
}
