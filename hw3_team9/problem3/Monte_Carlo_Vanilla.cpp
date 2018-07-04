/* Monte_Carlo_Vanilla.cpp
 * 09/17/17
 *
 * Chapin Day
 *
 * Uses Monte Carlo simulation to price a non-path-dependent European option
 *
 *
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <tuple>
#include "Monte_Carlo_Vanilla.hpp"

double Monte_Carlo_Vanilla::Price (Option& o, double sigma, double r, int num_steps, int p2) 
{
	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();
	int N = num_steps;

	double price{0};

	Random_Number_Generator RNG;
	std::shared_ptr<std::vector<double> > Z;
	//auto Z = RNG.Inverse_Transform(N); // generate normal vector
	switch(p2)
	{
		case 0:
			Z = RNG.Inverse_Transform(N);
			break;
		case 1:
			Z = RNG.Accept_Reject(N);
			break;
		case 2:
			Z = RNG.Box_Muller(N);
			break;
	}

	N = Z->size(); // reset # denominator to number of random numbers (differs for A-R & B-M)

	double Si, Vi;

	/*
	for (auto & z : *Z)
	{
		Si = S * pow(E,(r-q-sigma*sigma/2)*T + sigma*sqrt(T)*z); // UA price for path
		Vi = o.Payoff(Si); // option payoff for path
		//price += 1.0/N * Vi; // accumulate average option price
		price += Vi; // accumulate raw option price
		//std::cout << Si << '\t' << Vi << '\t' << price << std::endl;
	}
	*/
	auto iter1 = (*Z).begin();
	auto iter2 = (*Z).end();
	price = raw_Price(o, sigma, r, iter1, iter2);

	return price/N;
}

template <typename U>
double raw_Price(Option& o, double sigma, double r, U begin_iter, U end_iter)
{ // useful for multithreading
	// return raw sum of option prices calculated from a portion of the normal vector
	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();
	double Si, Vi, z;

	double price{0};

	int threshold = 500000;

	auto vector_length = end_iter - begin_iter;
	if (vector_length < threshold) // small enough to compute in a single thread
	{
		for (auto iter = begin_iter; iter != end_iter; ++iter)
		{ // iterate through portion of normal vector
			z = *iter;
			Si = S * pow(E,(r-q-sigma*sigma/2)*T + sigma*sqrt(T)*z); // UA price for path
			Vi = pow(E,-r*T) * o.Payoff(Si); // option payoff for path
			price += Vi; // accumulate raw option price
		}

		return price; // raw, unaveraged price
	}
	else
	{ // multithreaded approach for large vector
		auto new_begin = begin_iter + vector_length/2; // split vector in half
		std::future<double> fut = std::async(std::launch::async, [&]()
				{ return raw_Price(o,sigma,r,new_begin,end_iter); } );
		return raw_Price(o,sigma,r,begin_iter,new_begin) + fut.get();
	}
}	

std::tuple<double,double,double,double,double> Monte_Carlo_Vanilla::Greeks
(Option& o, double sigma, double r, int num_steps, int p2) 
{
	// delta,gamma,theta,vega, price

	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();
	int N = num_steps;

	double delta{0}, gamma{0}, theta{0}, vega{0}, price{0};

	Random_Number_Generator RNG;
	auto Z = RNG.Inverse_Transform(N); // generate normal vector

	auto iter1 = (*Z).begin();
	auto iter2 = (*Z).end();
	auto greeks = raw_Greeks(o, sigma, r, iter1, iter2);

	std::tie(delta,gamma,theta,vega,price) = greeks;

	return std::make_tuple(delta/N,gamma/N,theta/N,vega/N,price/N);
}

	template <typename U>
otuple raw_Greeks(Option& o, double sigma, double r, U begin_iter, U end_iter)
{ // useful for multithreading
	// return raw sum of option prices calculated from a portion of the normal vector
	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();
	double Si, Vi, z;

	double delta{0}, gamma{0}, theta{0}, vega{0}, price{0};

	/*
	// ERROR CHECKING
	double delta_i,vega_i,price_i,stock_i; 
	std::vector<std::tuple<double,double,double,double,double> > error_checking;
	*/

	int threshold = 100000;

	auto vector_length = end_iter - begin_iter;

	if (vector_length < threshold) // small enough to compute in a single thread
	{
		for (auto iter = begin_iter; iter != end_iter; ++iter)
		{ // iterate through portion of normal vector
			z = *iter;
			Si = S * pow(E,(r-q-sigma*sigma/2)*T + sigma*sqrt(T)*z); // UA price for path
			Vi = pow(E,-r*T)*o.Payoff(Si); // option payoff for path
			price += Vi; // accumulate raw option price
			delta += 0;
			vega += 0;
			//price_i=Vi; delta_i=0; vega_i=0; stock_i=Si; // ERROR CHECKING
			// indicator function: calculate Greeks for path if payoff is positive
			if (Vi > 0)
			{
				// recall o.type() is 1 for call, -1 for put
				delta += o.type() * pow(E,-r*T) * Si / S;
				vega += Si * o.type() * pow(E,-r*T) * (-sigma*T + sqrt(T) * z);
				//delta_i = o.type() * pow(E,-r*T) * Si / S;
				//vega_i = Si * o.type() * pow(E,-r*T) * (-sigma*T + sqrt(T) * z);
			}
			//error_checking.push_back(std::make_tuple(z,stock_i,price_i,delta_i,vega_i));
		}

		/*
		// ERROR CHECKING
		if (o.type() == 1)
		{
			std::ofstream myfile;
			myfile.open("vega_check.csv");
			for (auto & elem : error_checking)
			{
				std::tie(z,stock_i,price_i,delta_i,vega_i)=elem;
				myfile << z << ";" << stock_i << ";" << price_i << ";" << delta_i << ";" << vega_i << "\n";
			}
			myfile.close();
		}
		*/

		return std::make_tuple(delta,gamma,theta,vega,price); // raw, unaveraged price
	}
	else
	{ // multithreaded approach for large vector
		auto new_begin = begin_iter + vector_length/2; // split vector in half
		//std::cout << "multithreading " << vector_length/2 << " paths\n";
		std::future<otuple> fut = std::async(std::launch::async, [&]()
				{ return raw_Greeks(o,sigma,r,new_begin,end_iter); } );
		return tuple_add(raw_Greeks(o,sigma,r,begin_iter,new_begin), fut.get());
	}
}	

otuple tuple_add(const otuple& tup1, const otuple& tup2)
{ // add elements of two option greek result tuples
	double d1,d2,g1,g2,t1,t2,v1,v2,p1,p2;
	std::tie(d1,g1,t1,v1,p1) = tup1;
	std::tie(d2,g2,t2,v2,p2) = tup2;
	d1 += d2;
	g1 += g2;
	t1 += t2;
	v1 += v2;
	p1 += p2;
	return std::make_tuple(d1,g1,t1,v1,p1);
}
