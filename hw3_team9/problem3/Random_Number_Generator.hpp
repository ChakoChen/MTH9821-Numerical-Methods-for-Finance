/* Random_Number_Generator.hpp
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
 */

#ifndef RANDOM_NUMBER_GENERATOR_HPP
#define RANDOM_NUMBER_GENERATOR_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <tuple>

using pVec = std::shared_ptr<std::vector<double> >; // pointer to vector = return type

class Random_Number_Generator
{
	private:
		//pVec _random_nums;

	public:
		//Random_Number_Generator();							// Default ctor
		//Random_Number_Generator(const Random_Number_Generator&);				// Copy ctor
		//Random_Number_Generator(std::size_t); // Init ctor with length of vector
		//virtual ~Random_Number_Generator();					// Dtor
		//Random_Number_Generator& operator= (const Random_Number_Generator&);	// Assignment op

		pVec Uniform(std::size_t, unsigned long long=1,
				int=0, int=1); // return uniformly-distributed [0,1] (default if range unspecified)
		pVec Inverse_Transform(std::size_t, unsigned long long=1); // return standard normal 
		pVec Accept_Reject(std::size_t, unsigned long long=1); // return standard normal 
		pVec Box_Muller(std::size_t, unsigned long long=1);
};

double maximum(pVec p);
double minimum(pVec p);
double variance(pVec p);
double average(pVec p);
void stats(pVec p);



#endif // RANDOM_NUMBER_GENERATOR_HPP
