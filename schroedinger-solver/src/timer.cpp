#include "timer.hpp"

namespace solver {

Timer::Timer() : start(std::chrono::high_resolution_clock::now())
{
}

Timer::~Timer()
{
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  std::cout << "Execution time: " << duration << " μs\n";
}

};  // namespace solver