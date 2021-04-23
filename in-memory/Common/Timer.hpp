#ifndef __TIMER_HPP_
#define __TIMER_HPP_

#include <chrono>

using namespace std::chrono;
struct HighResolutionTimer {
  void restart() {
    _start = high_resolution_clock::now();
  }

  double elapsed() const {
    auto cur = high_resolution_clock::now();
    return duration<double, std::micro>(cur - _start).count();
  }

private:
  high_resolution_clock::time_point _start;
};

#endif // __TIMER_HPP_