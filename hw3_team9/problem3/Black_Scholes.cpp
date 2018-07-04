/* Black_Scholes.cpp
 * 09/04/17
 *
 * Chapin Day
 *
 * A child of the Option_Pricer class that employs the Black-Scholes method
 *
 */

#include <iostream>
#include "Black_Scholes.hpp"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.141592653589793238
#define E 2.718281828459045235
#define TOL 1e-12

using namespace std;


double i_func(double input_1)
{
    return 1 / sqrt(2 * PI) * pow(E, -0.5 * pow(input_1, 2));
}

// changheng's pdf of normal
double n(double x)
{
    return 1.0/sqrt(2.0 * M_PI) * exp(-x*x*0.5);
}


double integrate(double low_bound, double up_bound, int n)
{
    // simpsons rule
    double xi, ai;
    double step = (up_bound - low_bound)/static_cast<double>(n);
    double midpoint = step / 2;
    double result = 0;
    // term 1: h ( f(a0)/6 + f(an)/6)
    result += step / 6 * (i_func(up_bound) + i_func(low_bound));
    // term 2: h/3 * sum(i=1 -> n-1)f(ai)
    for (int i = 1; i < n; i++) {
	ai = low_bound + static_cast<double>(i) * step;	
	result += i_func(ai) * step / 3;
    }
    // term 3: 2h/3 * sum(i=1 -> n)f(xi)
    for (int i = 0; i < n; i++) {
	xi = (low_bound + static_cast<double>(i) * step + midpoint);
	result += i_func(xi) * 2 * step / 3;
    }
    return result;
}

double normal(double b, double tol)
{
    // b = upper bound
    double a = 0; // lower bound
    int n = 4; // intervals
    //double tol = 1e-12; // tolerance
    double result2, diff, diff2; 

    double result1 = 0.5 + integrate(a, b, n);

    //cout << "n\tintegral\tdiff" << endl;
    //cout << n << "\t" << result1 << endl;

    do {
	n = n * 2;
	result2 = 0.5 + integrate(a, b, n);
	diff = (result2 - result1);
	diff2 = (result1 - result2);
	if (diff2 > diff) diff = diff2;
	//cout << n << "\t" << result2 << "\t" << diff << endl;
	result1 = result2;
    } while (diff > tol);

    return result1;
}

// Approximation to the CDF from Changheng
double approx_normal(double x,double throw_away)
{
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}



/*
int main()
{
    double PCP = C_t - P_t - S * pow(E, -1 * q * (T-t))
	+ K * pow(E, -1 * r * (T-t));

    cout << "Put call parity: " << endl << "C_t - P_t: " << C_t - P_t 
	<< endl << "Se-qt - Ke-rt: " << S * pow(E, -1 * q * (T-t))
	- K * pow(E, -1 * r * (T-t)) << endl << "= " << PCP << endl << endl;

    return 0;
}
*/

double Black_Scholes::Price (Option& o, double sigma, double r, int num_steps, int p2) 
{ // num_steps only there to match function signature of Option_Pricer
	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();
	// r and sigma received...
	double result=0;

	double d1 = (log(S / K) + (r - q + pow(sigma, 2) / 2) * (T - t) ) / 
		(sigma * sqrt(T - t));
	double d2 = d1 - sigma * sqrt(T - t);

	if (o.type() == 1) // call
	{

		double Nd1 = normal(d1, TOL);
		double Nd2 = normal(d2, TOL);

		double cterm1 = S * Nd1 * pow(E, -1 * q * (T - t));
		double cterm2 = K * Nd2 * pow(E, -1 * r * (T - t));

		result = cterm1 - cterm2;
	}

    // put calculations

	else
	{ // must be a put

		double pNd2 = normal(-1 * d2, TOL);
		double pNd1 = normal(-1 * d1, TOL);
		double pterm1 = K * pow(E, -1 * r * (T - t)) * pNd2; 
		double pterm2 = S * pow(E, -1 * q * (T - t)) * pNd1; 

		result = pterm1 - pterm2;
	}

	return result;
}

std::tuple<double,double,double,double,double> Black_Scholes::Greeks (Option& o, double sigma, double r, int num_steps, int p2) 
{ // delta,gamma,theta,vega,price
    // call greeks
	double delta, gamma, theta, vega, price;

	double S = o.S();
	double K = o.K();
	double T = o.T();
	double t = 0;
	double q = o.q();

	double d1 = (log(S / K) + (r - q + pow(sigma, 2) / 2) * (T - t) ) / 
		(sigma * sqrt(T - t));
	double d2 = d1 - sigma * sqrt(T - t);


	if (o.type() == 1) // call
	{
		double Nd1 = normal(d1, TOL);
		double Nd2 = normal(d2, TOL);
		delta = pow(E, -1 * q * (T - t)) * Nd1;
		gamma = pow(E, -1 * q * (T - t)) 
			* pow(E, -1 * pow(d1, 2)/2)
			/ (S * sigma * sqrt(2 * PI * (T - t)));
		vega = S * pow(E, -1 * q * (T - t)) * sqrt(T - t)
			* pow(E, -1 * pow(d1, 2)/2)
			/ sqrt(2 * PI);
		theta = S * sigma * pow(E, -1 * q * (T - t))
			* pow(E, -1 * pow(d1, 2)/2)
			/ ( 2 * sqrt(2 * PI * (T - t)))
			+ q * S * pow(E, -1 * q * (T - t)) * Nd1
			- r * K * pow(E, -1 * r * (T - t)) * Nd2; 

		double cterm1 = S * Nd1 * pow(E, -1 * q * (T - t));
		double cterm2 = K * Nd2 * pow(E, -1 * r * (T - t));

		price = cterm1 - cterm2;
		//double c_rho = K * (T - t) * pow(E, -1 * r * (T - t)) * Nd2; 
	}
	else
	{ // must be a put

		// put greeks
		double pNd2 = normal(-1 * d2, TOL);
		double pNd1 = normal(-1 * d1, TOL);
		delta = -pow(E, -1 * q * (T - t)) * pNd1;
		gamma = pow(E, -1 * q * (T - t)) 
			* pow(E, -1 * pow(d1, 2)/2)
			/ (S * sigma * sqrt(2 * PI * (T - t)));
		vega = S * pow(E, -1 * q * (T - t)) * sqrt(T - t)
			* pow(E, -1 * pow(d1, 2)/2)
			/ sqrt(2 * PI);
		theta = -S * sigma * pow(E, -1 * q * (T - t))
			* pow(E, -1 * pow(d1, 2)/2)
			/ ( 2 * sqrt(2 * PI * (T - t)))
			- q * S * pow(E, -1 * q * (T - t)) * pNd1
			+ r * K * pow(E, -1 * r * (T - t)) * pNd2; 



		double pterm1 = K * pow(E, -1 * r * (T - t)) * pNd2; 
		double pterm2 = S * pow(E, -1 * q * (T - t)) * pNd1; 

		price = pterm1 - pterm2;
		//double p_rho = -1 * K * (T - t) * pow(E, -1 * r * (T - t)) * pNd2; 
	}

	return std::make_tuple(delta, gamma, theta, vega, price);
}


