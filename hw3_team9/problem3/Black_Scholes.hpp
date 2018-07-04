/* Black_Scholes.hpp
 * 09/04/17
 *
 * Chapin Day
 *
 * An Option_Pricer class that uses the Black-Scholes method
 *
 */

#ifndef BLACK_SCHOLES_HPP
#define BLACK_SCHOLES_HPP

#include <iostream>
#include "Option.hpp"
#include "Option_Pricer.hpp"

class Black_Scholes : public Option_Pricer
{
	private:

	public:
		/*
		Black_Scholes() ;							// Default ctor
		Black_Scholes(const Black_Scholes&);				// Copy ctor
		Black_Scholes& operator= (const Black_Scholes&);	// Assignment op
		virtual ~Black_Scholes() ;					// Dtor
		*/
	
		virtual double Price (Option& o, double sigma, double r, int num_steps=0,int=0) ;
		virtual std::tuple<double,double,double,double,double> Greeks
			(Option& o, double sigma, double r, int num_steps=0,int=0) ; // delta,gamma,theta,vega
};

#endif // BLACK_SCHOLES_HPP
