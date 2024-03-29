/**
 * @file fluctuation.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef INITIALIZATION_FLUCTUATION_HPP_
#define INITIALIZATION_FLUCTUATION_HPP_

#include <omp.h>

#include <array>
#include <iostream>
#include <random>
#include <vector>

#include "brown_factor.hpp"

class fluctuation : virtual public brown_factor {
 private:
  /* data */

 protected:
  std::vector<std::array<double, 3>> g0;
  std::vector<std::array<double, 3>> g1;

 public:
  void init_fluctuation();
  void generate_Gamma();
};

#endif  // INITIALIZATION_FLUCTUATION_HPP_
