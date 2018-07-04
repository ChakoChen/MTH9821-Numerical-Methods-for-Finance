/* Monte_Carlo_Path_Dependent.hpp
 * 09/17/17
 *
 * Chapin Day
 *
 * This will return the value of an option that depends on a lognormally-distributed
 * underlying asset.  *
 */

#ifndef MONTE_CARLO_PATH_DEPENDENT_HPP
#define MONTE_CARLO_PATH_DEPENDENT_HPP

#include <iostream>
#include "Random_Number_Generator.hpp"
#include "Option_Pricer.hpp"
#include "Option.hpp"
#include <memory>
#include <functional>
#include <cmath>
#include <tuple>
#include <thread>
#include <mutex>
#include <future>
#define E 2.718281828459045235


class Monte_Carlo_Path_Dependent : public Option_Pricer
{
	private:

	public:
		/*
		Monte_Carlo_Path_Dependent();							// Default ctor
		Monte_Carlo_Path_Dependent(const Monte_Carlo_Path_Dependent&);				// Copy ctor
		virtual ~Monte_Carlo_Path_Dependent();					// Dtor
		Monte_Carlo_Path_Dependent& operator= (const Monte_Carlo_Path_Dependent&);	// Assignment op
		*/
		double Price (Option& o, double sigma, double r, int num_steps, int param2=200) ;
		std::tuple<double,double,double,double,double> Greeks
			(Option& o, double sigma, double r, int num_steps, int param2=200) ; // delta,gamma,theta,vega

};

#endif // MONTE_CARLO_PATH_DEPENDENT_HPP
