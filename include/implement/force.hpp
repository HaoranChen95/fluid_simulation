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

#ifndef FORCE_HPP_
#define FORCE_HPP_

#include "position.hpp"

class force : protected position {
 private:
  double E_pot_;
  double LJ(uint64_t i, uint64_t j);

 protected:
  std::vector<std::array<double, 3>> f0;
  std::vector<std::array<double, 3>> f1;

 public:
  void init_force();
  void calc_force();
  double E_pot() const;
  ~force();
};

#endif  // FORCE_HPP_