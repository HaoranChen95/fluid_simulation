/**
 * @file position.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef IMPLEMENTATION_POSITION_HPP_
#define IMPLEMENTATION_POSITION_HPP_

#include <omp.h>

#include <array>
#include <iostream>
#include <random>
#include <vector>

#include "box.hpp"

class position : virtual public box {
 private:
 protected:
  std::vector<std::array<double, 3>> r;
  std::vector<std::array<double, 3>> dr;

 public:
  void init_position(/* args */);
  double minium_image(const uint64_t &i, const uint64_t &j, const int &ax);
  ~position();
};

#endif  // IMPLEMENTATION_POSITION_HPP_
