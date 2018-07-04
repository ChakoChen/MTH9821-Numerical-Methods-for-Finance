/* Option.hpp
 * 09/04/17
 *
 * Chapin Day
 *
 * Note: option arguments: American(T/F), type (1/-1), S, K, T, q, Barrier, barrier_type
 * where barrier_type is 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in
 *
 */

#ifndef OPTION_HPP
#define OPTION_HPP

#include <iostream>

class Option
{
	private:
		double _K; // strike
		double _S; // underlying price
		double _T; // maturity
		double _q; // dividend rate
		int _type; // 1 = call; -1 = put
		bool _American;
		double _B; // barrier
		int _barrier_type; // 1=up-and-out, 2=up-and-in, 3=down-and-out, 4=down-and-in

	public:
		Option();							// Default ctor
		Option(bool Amer, int type, double s, double k, double t, double q, double B=0, int bt=0)
			: _American(Amer), _K(k), _S(s), _T(t), _q(q), _type(type),
			_B(B), _barrier_type(bt) {} // init ctor
		Option(const Option&);				// Copy ctor
		Option& operator= (const Option&);	// Assignment op
		virtual ~Option() {};					// Dtor

		// getters, setters
		int type() { return _type; }
		double S() { return _S; }
		double K() { return _K; }
		double T() { return _T; }
		double q() { return _q; }
		double B() { return _B; }
		int barrier_type() { return _barrier_type; }
		bool American() { return _American; }
		void S(double s) { _S = s; }
		void B(double b) { _B = b; }
		void American_Toggle() { _American = !_American; } // toggle American
		double Payoff(double S);

		friend std::ostream& operator << (std::ostream&, const Option&);
};


#endif // OPTION_HPP
