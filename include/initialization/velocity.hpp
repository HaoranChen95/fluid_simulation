/**
 * @file velocity.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef INITIALIZATION_VELOCITY_HPP_
#define INITIALIZATION_VELOCITY_HPP_

#include <omp.h>

#include <array>
#include <iostream>
#include <random>
#include <vector>

#include "box.hpp"

class velocity : virtual public box {
 private:
  double E_kin_;
  double mean_E_kin_ = 0;
  uint64_t step_ = 0;
  void calc_mean_E_kin();

 protected:
  std::vector<std::array<double, 3>> v;

 public:
  void init_velocity();
  void vel_correcter();
  void calc_E_kin();
  double E_kin() const;
  ~velocity();
};

#endif  // INITIALIZATION_VELOCITY_HPP_
