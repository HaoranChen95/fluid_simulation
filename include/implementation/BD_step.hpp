/**
 * @file BD_step.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief 
 * @version 0.1
 * @date 2021-12-26
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef IMPLEMENTATION_BD_STEP_HPP_
#define IMPLEMENTATION_BD_STEP_HPP_

#include "initialization.hpp"

class BD_step : virtual protected initialization{
 private:
  /* data */
  void calc_vel();
  void calc_pos();

 public:
  BD_step();
  void run_BD_step();
  ~BD_step();
};

#endif  // IMPLEMENTATION_BD_STEP_HPP_
