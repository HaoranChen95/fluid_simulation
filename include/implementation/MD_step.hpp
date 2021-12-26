/**
 * @file MD_step.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef IMPLEMENTATION_MD_STEP_HPP_
#define IMPLEMENTATION_MD_STEP_HPP_

#include "initialization.hpp"

class MD_step : virtual protected initialization{
 private:
  /* data */
  void calc_vel();
  void calc_pos();

 public:
  MD_step();
  void run_MD_step();
  ~MD_step();
};

#endif  // IMPLEMENTATION_MD_STEP_HPP_
