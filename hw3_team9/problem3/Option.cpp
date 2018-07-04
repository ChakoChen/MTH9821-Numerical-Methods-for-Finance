/* Option.cpp
 * 09/04/17
 *
 * Chapin Day
 *
 * just to contain the print function?
 *
 */

#include <iostream>
#include "Option.hpp"

Option::Option(const Option& source)
{
	_K = source._K; // strike
	_S = source._S; // underlying price
	_T = source._T; // maturity
	_q = source._q; // dividend rate
	_type = source._type; // 1 = call; -1 = put
	_American = source._American;
	_B = source._B; // barrier
	_barrier_type = source._barrier_type; // 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in
}


std::ostream& operator << (std::ostream& os, const Option& o)
{
	if (o._American) os << "American ";
	if (o._B > 0)
	{
		// 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in
		switch(o._barrier_type)
		{
			case 1:
				os << "Up-and-Out "; break;
			case 2:
				os << "Up-and-In "; break;
			case 3:
				os << "Down-and-Out "; break;
			case 4:
				os << "Down-and-In "; break;
			default:
				os << "ERROR - no barrier definition ";
		}
	}

	if (o._type == 1 ) os << "Call: ";
	else os << "Put: ";
	os << "K: " << o._K << ", S: " << o._S << ", T: " << o._T << ", q: " << o._q;
	if (o._B > 0) os << ", B: " << o._B;
	return os;
}	

double Option::Payoff(double S)
{
	double payoff = _type * (S - _K ); // S-K for call, K-S for put
	return ( payoff > 0 ? payoff : 0 );
}
