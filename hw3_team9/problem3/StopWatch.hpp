/* StopWatch.hpp
 * 7/14/17
 *
 * An adapter class for certain functionality in std::chrono.
 * 
 * It times how long it takes for an operation to complete.
 */

#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include <chrono>
#include <ratio>

class StopWatch
{
	public:
		StopWatch(); // def ctor

		void StartStopWatch();
		void StopStopWatch();
		void Reset();

		double GetTime() const;
		int GetClicks() const;

	private:
		StopWatch (const StopWatch&); // copy ctor "deleted"
		StopWatch& operator= (const StopWatch&); // assignment op "deleted"
		std::chrono::steady_clock::time_point _start; // starting time point
		std::chrono::duration<unsigned long long,std::nano> _dur; // duration of event
		//decltype(std::chrono::steady_clock::period::den) _den; // denom of time period
		unsigned long long _den; // denom of time period
};


#endif // STOPWATCH_HPP
