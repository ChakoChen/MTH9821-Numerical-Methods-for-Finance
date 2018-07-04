/* Monte_Carlo_Vanilla.hpp
 * 09/17/17
 *
 * Chapin Day
 *
 * Uses Monte Carlo simulation to price a non-path-dependent European option
 *
 * param2 may be used to select random number generation method:
 * 0 = Inverse Transform
 * 1 = Acceptance-Rejection
 * 2 = Box-Muller
 *
 */

#ifndef MONTE_CARLO_VANILLA_HPP
#define MONTE_CARLO_VANILLA_HPP

#include <iostream>
#include "Random_Number_Generator.hpp"
#include "Option_Pricer.hpp"
#include "Option.hpp"
#include <cmath>
#include <tuple>
#include <thread>
#include <mutex>
#include <future>
#define E 2.718281828459045235

using otuple = std::tuple<double,double,double,double,double>;


class Monte_Carlo_Vanilla : public Option_Pricer
{
	private:

	public:
		/*
		Monte_Carlo_Vanilla();							// Default ctor
		Monte_Carlo_Vanilla(const Monte_Carlo_Vanilla&);				// Copy ctor
		virtual ~Monte_Carlo_Vanilla();					// Dtor
		Monte_Carlo_Vanilla& operator= (const Monte_Carlo_Vanilla&);	// Assignment op
		*/
		double Price (Option& o, double sigma, double r, int num_steps, int param2=0) ;
		std::tuple<double,double,double,double,double> Greeks
			(Option& o, double sigma, double r, int num_steps, int =0) ; // delta,gamma,theta,vega
};

template <typename U>
double raw_Price(Option& o, double sigma, double r, U begin_iter, U end_iter);
template <typename U>
otuple raw_Greeks(Option& o, double sigma, double r, U begin_iter, U end_iter);
otuple tuple_add(const otuple&, const otuple&);

#endif // MONTE_CARLO_VANILLA_HPP
