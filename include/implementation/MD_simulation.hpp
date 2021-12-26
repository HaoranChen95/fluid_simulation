/**
 * @file MD_SIMULATION.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef IMPLEMENTATION_MD_SIMULATION_HPP_
#define IMPLEMENTATION_MD_SIMULATION_HPP_

#include "initialization.hpp"
#include "output.hpp"

class MD_simulation : virtual protected initialization, virtual public output {
 private:
  /* data */
  void calc_vel();
  void calc_pos();

 public:
  MD_simulation();
  void MD_relaxation();
  void MD_implementation();
  ~MD_simulation();
};

#endif  // IMPLEMENTATION_MD_SIMULATION_HPP_
