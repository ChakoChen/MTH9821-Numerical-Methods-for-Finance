/* Option_Pricer.hpp
 * 09/04/17
 *
 * Chapin Day
 *
 * Base class for option pricers
 *
 */

#ifndef OPTION_PRICER_HPP
#define OPTION_PRICER_HPP

#include <iostream>
#include "Option.hpp"
#include <tuple>

class Option_Pricer
{
	private:

	public:
		/*
		Option_Pricer();							// Default ctor
		Option_Pricer(const Option_Pricer&);				// Copy ctor
		Option_Pricer& operator= (const Option_Pricer&);	// Assignment op
		virtual ~Option_Pricer();					// Dtor
		*/

		virtual double Price (Option& o, double sigma, double r, int num_steps=0, int param2=0) =0;
		virtual std::tuple<double,double,double,double,double> Greeks
			(Option& o, double sigma, double r, int num_steps=0, int param2=0) =0; // delta,gamma,theta,vega
		std::tuple<double,int> Tolerance_Price (Option& o, double sigma, double r, int num_steps,
				double tolerance)
		{
			double diff, px_new;
			std::cout << "Tolerance pricing to tolerance: " << tolerance << std::endl;
			double px_old = Price(o, sigma, r, num_steps,0);
			do
			{
				num_steps *= 2;
				px_new = Price(o, sigma, r, num_steps,0);
				//std::cout << num_steps << " " << px_new << std::endl;
				diff = px_new - px_old;
				diff = ( diff < 0 ? -diff : diff );
				px_old = px_new;
			} while (diff > tolerance);
			return std::make_tuple(px_new,num_steps);
		}

};

#endif // OPTION_PRICER_HPP
