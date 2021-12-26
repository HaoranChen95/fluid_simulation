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

#ifndef IMPLEMENTATION_VELOCITY_HPP_
#define IMPLEMENTATION_VELOCITY_HPP_

#include <omp.h>

#include <array>
#include <iostream>
#include <random>
#include <vector>

class velocity {
 private:
  /* data */
 protected:
  std::vector<std::array<double, 3>> v;

 public:
  void init_velocity(const uint64_t Nm, const double kT);
  ~velocity();
};

#endif  // IMPLEMENTATION_VELOCITY_HPP_
