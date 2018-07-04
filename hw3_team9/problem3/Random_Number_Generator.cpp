/* Random_Number_Generator.cpp
 * 09/17/17
 *
 * Chapin Day
 *
 * This generates a vector of uniform random numbers between zero and 1
 * using the Linear Congruential Generator, and then transforms it into 
 * a vector of standard-normally-distributed random numbers using either:
 * - Inverse Transform Method
 * - Acceptance-Rejection Method
 * - Box-Muller Method
 *
 * It returns a pointer to a vector of random numbers, of length specified by user
 *
 *
 */

#include "Random_Number_Generator.hpp"
#define E 2.718281828459045235

pVec Random_Number_Generator::Uniform(std::size_t length, unsigned long long seed,
		int lower_bound, int upper_bound)
{
	// return uniformly-distributed [0,1] random numbers
	pVec result (new std::vector<double>); // output, pointer to vector of doubles
	
	// Linear Congruential Generator:
	// Choose positive integers: X0, a, c, m 
	// X_i+1 = (a*x_i + c)mod m
	// U_i+1 = X_i+1 / m
	double range = 1.0 * (upper_bound - lower_bound); // defaults to 1.0
	
	unsigned long long X0{seed}, a{39373}, c{0}, m(pow(2,31)-1);

	unsigned long long X_old=X0;

	for (int i=0; i<length; ++i)
	{
		/*
		std::cout << "num: " << a << " * " << X_old << " + " << c << " = " 
			<< a * X_old + c << "; den: " << m << std::endl;
			*/
		unsigned long long X_new = (a * X_old + c) % m;
		result->push_back( (range / m) * X_new  + lower_bound);
		X_old = X_new;
	}


	return result;
}

pVec Random_Number_Generator::Inverse_Transform(std::size_t length, unsigned long long seed)
{
	// return standard normal random numbers
	pVec result = Uniform(length,seed); // first use LCG to generate Uniform
	// constants from website: http://www.unc.edu/~gfish/cndev.c
	double A[4]={
		2.50662823884,
		-18.61500062529,
		41.39119773534,
		-25.44106049637
	};
	double B[4]={
		-8.47351093090,
		23.08336743743,
		-21.06224101826,
		3.13082909833
	};
	double C[9]={
		0.3374754822726147,
		0.9761690190917186,
		0.1607979714918209,
		0.0276438810333863,
		0.0038405729373609,
		0.0003951896511919,
		0.0000321767881768,
		0.0000002888167364,
		0.0000003960315187
	};
	// Method from Glasserman, p. 68, fig. 2.13, approximates inverse normal:
	// aka Beasley-Springer-Moro algorithm 
	double x,y,r;
	for (auto & u : *result) // transform uniform vector, element by element
	{
		y = u - 0.5;
		if (std::abs(y) < 0.42)
		{
			r = y * y;
			x = y * (((A[3]*r + A[2])*r + A[1])*r + A[0]) /
				((((B[3]*r + B[2])*r + B[1])*r + B[0])*r + 1);
		}
		else
		{
			r = u;
			if (y > 0) r = 1-u;
			r = log(-log(r));
			x = C[0]+r*(C[1] + r*(C[2] + r*(C[3] +r*(C[4] 
								+ r*(C[5] + r*(C[6] +r*(C[7] + r*C[8])))))));
			if (y < 0) x = -x;
		}
		u = x; // replace uniform with normal 
	}
	return result;
}


pVec Random_Number_Generator::Accept_Reject(std::size_t length, unsigned long long seed)
{
	// return standard normal random numbers

	// need 3 uniforms for each normal
	pVec result (new std::vector<double> ); // output vector
	double x, u1, u2, u3;

	/* creates "length" normals; homework says this is wrong
	do
	{
		// create 3 uniform random variables [0,1]
		pVec U = Uniform(3,seed);
		auto u_iter = (*U).begin();
		u1 = *u_iter++;
		u2 = *u_iter++;
		u3 = *u_iter;

		x = -log(u1);
		if (u2 > pow(E,-0.5*(x-1)*(x-1))) 
		{
			seed = u3 * m;	// reset seed so don't repeat same three random variables
			//std::cout << seed << std::endl;
			continue; // go back to start of loop
		}
		if (u3 <= 0.5) x = -x; // reverse sign randomly
		result->push_back(x); // enter result into normal vector
		seed = u3 * m;	// reset seed so don't repeat same three random variables
		++count; // increment result count
	} while (count < length);
	*/

	pVec U = Uniform(length,seed); // N uniforms; we'll print out how many normals are returned

	int max_normals = length/3; // must work from independent groups of 3 uniforms
	//std::cout << "maximum normals: " << max_normals << std::endl;
	auto u_iter = (*U).begin();

	for (int i=0; i<max_normals; ++i)
	{
		u1 = *u_iter++;
		u2 = *u_iter++;
		u3 = *u_iter++;
		x = -log(u1);
		if (u2 > pow(E,-0.5*(x-1)*(x-1))) 
		{
			continue; // reject!
		}
		if (u3 <= 0.5) x = -x; // reverse sign randomly
		result->push_back(x); // enter result into normal vector
	}

	std::cout << result->size() << " normals;";

	return result;
}

pVec Random_Number_Generator::Box_Muller(std::size_t length, unsigned long long seed)
{
	// return standard normal random numbers

	// need 2 uniforms for each normal
	pVec result (new std::vector<double> ); // output vector
	double x, y, u1, u2, z1, z2;

	/* old way, generated length normals; homework says this is wrong
	do 
	{
		pVec U = Uniform(2,seed,-1,1); // 2 uniform[-1,1] numbers
		auto u_iter = (*U).begin();
		u1=*u_iter++;
		u2=*u_iter;

		x = u1*u1 + u2*u2;
		if (x >= 1.0) // point is outside unit circle, so reject
		{
			seed = u2 * m; // set new seed
			continue;
		}

		y = sqrt(-2.0*log(x)/x);
		z1 = u1*y;
		z2 = u2*y;
		result->push_back(z1);
		result->push_back(z2);
		seed = u2 * m; // set new seed
		++count; ++count;
	} while (count < length);

	if (count > length) // added 2 z's per iteration; could have overshot
	{
		result->resize(length);
	}
	*/

	pVec U = Uniform(length,seed,-1,1); // generate N uniforms
	auto u_iter = (*U).begin();
	int max_length = length/2; // maximum # of b-m normals; must be even
	for (int i=0; i<max_length; ++i)
	{
		u1 = *u_iter++;
		u2 = *u_iter++;
		x = u1*u1 + u2*u2;
		if (x >= 1.0) // point is outside unit circle, so reject
		{
			continue;
		}

		y = sqrt(-2.0*log(x)/x);
		z1 = u1*y;
		z2 = u2*y;
		result->push_back(z1);
		result->push_back(z2);
	}

	std::cout << result->size() << " normals;";

	return result;

}

///////////////////////////////////// Vector Statistics /////////////////
double maximum(pVec p)
{ // return the maximum element of a vector
	double max = 0;
	for (auto & elem : *p)
	{
		if (elem > max) max = elem;
	}
	return max;
}

double minimum(pVec p)
{ // return the minimum element of a vector
	double min = 0;
	for (auto & elem : *p)
	{
		if (elem < min) min = elem;
	}
	return min;
}

double average(pVec p)
{ // return the average element of a vector
	double avg = 0;
	double sum = 0;
	int count = 0;
	for (auto & elem : *p)
	{
		++count;
		sum += elem;
	}
	avg = sum / count;
	return avg;
}

double variance(pVec p)
{
	double mean = average(p);
	double variance = 0;
	for (auto & elem : *p)
	{
		variance += pow(elem - mean, 2);
	}
	variance /= (*p).size();
	return variance;
}

void stats(pVec p)
{
	std::cout << "Length: " << (*p).size()
		<< "\nMin: " << minimum(p)
		<< "\nMax: " << maximum(p)
		<< "\nAvg: " << average(p)
		<< "\nVar: " << variance(p)
		<< std::endl;
}
