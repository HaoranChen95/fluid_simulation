/**
 * @file fluid_simulation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "fluid_simulation.hpp"

fluid_simulation::fluid_simulation(const int argc, const char **argv) {
  init(argc, argv);
}

void fluid_simulation::relax() {
  if (gamma()) {
    BD_relaxation();
  } else {
    MD_relaxation();
  }
  std::cout << "finished relation" << std::endl;
}
void fluid_simulation::implement() {
  if (gamma()) {
    BD_implementation();
  } else {
    MD_implementation();
  }
  std::cout << "finished implementation" << std::endl;
}

fluid_simulation::~fluid_simulation() {}
