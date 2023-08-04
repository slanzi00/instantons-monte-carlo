#ifndef TIMER_HPP
#define TIMER_HPP
#define ENABLE_TIME_MEASUREMENT

#include <chrono>
#include <iostream>

namespace solver {

class Timer
{
  std::chrono::time_point<std::chrono::high_resolution_clock> start;

 public:
  Timer();
  ~Timer();
};

};  // namespace solver

#endif