/**
 * @file brown_factor.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "brown_factor.hpp"

double brown_factor::C(const double x) {
  long double x_ = static_cast<long double>(x);
  return 2. * x_ - 3. + 4. * exp(-x_) - exp(-2. * x_);
}
double brown_factor::G(const double x) {
  long double x_ = static_cast<long double>(x);
  return exp(x_) - 2. * x_ - exp(-x_);
}

void brown_factor::calc_BD_factor() {
  double gh = gamma() * h();

  double C_gh = C(gh);
  double G_gh = G(gh);
  double E_gh = -C(gh) * C(-gh) - pow(G(gh), 2.);

  std::cout << "gh\t" << gh << "\tC\t" << C_gh << "\tG\t" << G_gh << "\tE\t"
            << E_gh << std::endl;
  BD_g0_1_ = sqrt(kT() / m() / gamma() / gamma() * C_gh);
  BD_g1_1_ = sqrt(kT() / m() / gamma() / gamma() * E_gh / C_gh);
  BD_g1_2_ = G_gh / C_gh;

  BD_r_1_ = (1. - exp(-gh)) / gamma();
  BD_r_2_ = (gh - 1. + exp(-gh)) / gamma() / gamma();

  BD_v_1_ = gamma() / (exp(gh) - 1.);
  BD_v_2_ = (gh - 1. + exp(-gh)) / gamma() / (exp(gh) - 1.);
  BD_v_3_ = (-gh - 1. + exp(gh)) / gamma() / (exp(gh) - 1.);
}

double brown_factor::BD_r_1() const { return BD_r_1_; }
double brown_factor::BD_r_2() const { return BD_r_2_; }
double brown_factor::BD_v_1() const { return BD_v_1_; }
double brown_factor::BD_v_2() const { return BD_v_2_; }
double brown_factor::BD_v_3() const { return BD_v_3_; }
double brown_factor::BD_g0_1() const { return BD_g0_1_; }
double brown_factor::BD_g1_1() const { return BD_g1_1_; }
double brown_factor::BD_g1_2() const { return BD_g1_2_; }

void brown_factor::print_brown_factor() {
  std::cout << "====== brown factor ======" << std::endl;
  std::cout << "BD_r_1\t" << BD_r_1_ << "\tBD_r_2\t" << BD_r_2_ << std::endl
            << "BD_v_1\t" << BD_v_1_ << "\tBD_v_2\t" << BD_v_2_ << "\tBD_v_3\t"
            << BD_v_3_ << std::endl
            << "BD_g0_1\t" << BD_g0_1_ << std::endl
            << "BD_g1_1\t" << BD_g1_1_ << "\tBD_g1_2\t" << BD_g1_2_
            << std::endl;
}
