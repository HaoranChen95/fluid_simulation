/**
 * @file box.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "box.hpp"

void box::kT(const double input) { kT_ = input; }
double box::kT() const { return kT_; }

/**
 * @brief set the periodic boundary condition
 *
 * @param l_b_x
 * @param l_b_y
 * @param l_b_z
 */
void box::l_b(const int ax, const double input) {
  l_b_[ax] = input;
  half_l_b_[ax] = 0.5 * input;
  inv_l_b_[ax] = 1. / input;
}
std::array<double, 3> box::l_b() const { return l_b_; }
std::array<double, 3> box::half_l_b() const { return half_l_b_; }
std::array<double, 3> box::inv_l_b() const { return inv_l_b_; }

void box::Nm(const uint64_t input) { Nm_ = input; }
void box::density(const double input) { density_ = input; }

void box::calc_Nm() {
  Nm_ = static_cast<uint64_t>(l_b_[0] * l_b_[1] * l_b_[2] * density_ / M_PI_4 /
                              sig2());
}

void box::calc_density() {
  density_ =
      static_cast<double>(Nm_) * M_PI_4 * sig2() / l_b_[0] / l_b_[1] / l_b_[2];
}
uint64_t box::Nm() const { return Nm_; }
double box::density() const { return density_; }
