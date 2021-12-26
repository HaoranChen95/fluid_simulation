/**
 * @file fluid_simulation.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef FLUID_SIMULATION_HPP_
#define FLUID_SIMULATION_HPP_

#include "implementation.hpp"
#include "initialization.hpp"
#include "output.hpp"

class fluid_simulation : virtual public initialization,
                         virtual public implementation,
                         virtual public output {
 private:
  /* data */
 public:
  fluid_simulation(const int argc, const char **argv);
  ~fluid_simulation();
};

#endif  // FLUID_SIMULATION_HPP_
