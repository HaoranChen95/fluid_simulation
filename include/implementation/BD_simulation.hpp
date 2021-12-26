/**
 * @file BD_SIMULATION.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef IMPLEMENTATION_BD_SIMULATION_HPP_
#define IMPLEMENTATION_BD_SIMULATION_HPP_

#include "initialization.hpp"
#include "output.hpp"

class BD_simulation : virtual protected initialization, virtual public output {
 private:
  /* data */
  void calc_vel();
  void calc_pos();

 public:
  BD_simulation();
  void BD_relaxation();
  void BD_implementation();
  ~BD_simulation();
};

#endif  // IMPLEMENTATION_BD_SIMULATION_HPP_
