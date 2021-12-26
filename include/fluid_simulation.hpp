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

#include "BD_simulation.hpp"
#include "MD_simulation.hpp"

class fluid_simulation : virtual public MD_simulation,
                         virtual public BD_simulation {
 private:
  /* data */
 public:
  fluid_simulation(const int argc, const char **argv);
  void relax();
  void implement();
  ~fluid_simulation();
};

#endif  // FLUID_SIMULATION_HPP_
