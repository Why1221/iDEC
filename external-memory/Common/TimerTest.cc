#include "Timer.hpp"
#include <iostream>
#include <thread>

int main() {

  using namespace std::chrono_literals;

  HighResolutionTimer timer;

  timer.restart();

  std::this_thread::sleep_for(3s);

  auto el = timer.elapsed();

  fprintf(stdout, "exec time: %.6f us\n", el);

  return 0;
}