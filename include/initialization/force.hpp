/**
 * @file force.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef INITIALIZATION_FORCE_HPP_
#define INITIALIZATION_FORCE_HPP_

#include <omp.h>

#include <array>
#include <iostream>
#include <random>
#include <vector>

#include "position.hpp"
#include "cell_list.hpp"

class force : virtual protected position , public cell_list {
 private:
  double E_pot_;
  void LJ(uint64_t i, uint64_t j);

 protected:
  std::vector<std::array<double, 3>> f0;
  std::vector<std::array<double, 3>> f1;

 public:
  void init_force();
  void calc_force();
  double E_pot() const;
};

#endif  // INITIALIZATION_FORCE_HPP_
