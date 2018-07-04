//This library contains Bisection/Newton's/Secant method to solve for x in f(x)=0


#ifndef EQUATIONSOLVER_HPP
#define EQUATIONSOLVER_HPP

#include <functional>
#include <cmath>

double Bisection(double a, double b, std::function<double(double)> fn)
{
	double tol_approx = 0.000001;
	double tol_consec = 0.000001;

	double xm;

	double xl = a, xr = b;
	double fl = fn(xl), fr = fn(xr);
	while ((std::abs(fl) > std::abs(fr) ? std::abs(fl) : std::abs(fr)) > tol_approx || (xr - xl) > tol_consec)
	{
		xm = (xl + xr) / 2;
		fl = fn(xl);
		fr = fn(xr);
		if (fl*fr < 0)
			xr = xm;
		else xl = xm;
	}
	return xm;
}

double Newton(double x0, std::function<double(double)> fn, std::function<double(double)> dfn)
{//Use Newton's method to solve f(x)=0
	double tol_approx = 0.000001;	//Newton method stopped at (1) |f(x_new)|<=tol_approx
	double tol_consec = 0.000001;	//Newton method stopped at (2) |x_new-x_old|<=tol_consec

									//set initial guesses
	double x_new = x0;
	double x_old = x0 - 1;

	//The recursion is x_k+1 = x_k - f(x_k)/f'(x_k)
	while (std::abs(fn(x_new)) > tol_approx || std::abs(x_new - x_old) > tol_consec)
	{
		x_old = x_new;
		x_new = x_old - fn(x_old) / dfn(x_old);
	}

	return x_new;
}



double Secant(double x0, double x1, std::function<double(double)> fn)
{//Use Secant's method to solve f(x)=0
	double tol_approx = 0.000001;	//Secant method stopped at (1) |f(x_new)|<=tol_approx
	double tol_consec = 0.000001;	//Secant method stopped at (1) |x_new-x_old|<=tol_consec

									//set initial guesses
	double x_new = x1;
	double x_old = x0;
	double x_oldest = 0;

	//The recursion is x_k+1 = x_k - (x_k-x_k-1)f(x_k)/(f(x_k)-f(x_k-1))
	while (std::abs(fn(x_new)) > tol_approx || std::abs(x_new - x_old) > tol_consec)
	{
		x_oldest = x_old;
		x_old = x_new;
		x_new = x_old - fn(x_old)*(x_old - x_oldest) / (fn(x_old) - fn(x_oldest));
	}

	return x_new;
}


#endif
