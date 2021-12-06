/**
 * @file box.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef BOX_HPP_
#define BOX_HPP_

#include <array>
#include <cstdint>

#include "particle_parameter.hpp"

class box : public particle_parameter {
 private:
  double kT_;
  double density_;
  std::array<double, 3> l_b_;
  std::array<double, 3> half_l_b_;
  std::array<double, 3> inv_l_b_;
  double Nm_;

 public:
  void kT(const double input);
  double kT() const;

  void l_b(const int ax, const double input);
  std::array<double, 3> l_b() const;
  std::array<double, 3> half_l_b() const;
  std::array<double, 3> inv_l_b() const;

  void Nm(const uint64_t input);
  void density(const double input);
  uint64_t Nm() const;
  void calc_Nm();
  void calc_density();
  double density() const;
};

#endif  // BOX_HPP_
