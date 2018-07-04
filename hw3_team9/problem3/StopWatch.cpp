/* StopWatch.cpp
 * 7/14/17
 *
 * Chapin Day
 *
 * An adapter class for certain functionality in std::chrono.
 * 
 * It times how long it takes for an operation to complete.
 */

#include "StopWatch.hpp"

#include <chrono>
#include <iostream>
#include <ratio>

StopWatch::StopWatch()
{ // def ctor
	_den = std::chrono::steady_clock::period::den;
}

void StopWatch::StartStopWatch()
{ // begin timing
	_start = std::chrono::steady_clock::now();
}

void StopWatch::StopStopWatch()
{ // end timing, calculate duration (clicks)
	_dur = std::chrono::steady_clock::now() - _start;
}

double StopWatch::GetTime() const
{ // return double of duration (milliseconds)
	unsigned long long numClicks = _dur.count();
	double seconds = double(numClicks) / _den;
	double milliseconds = seconds * 1'000;
	return milliseconds;
}

int StopWatch::GetClicks() const
{ // return system clock clicks
	return _dur.count();
}

void StopWatch::Reset()
{ // reset stopwatch; is this really necessary?
	_dur = std::chrono::duration<long,std::nano> ::zero();
}

