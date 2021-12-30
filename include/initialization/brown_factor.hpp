/**
 * @file brown_factor.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef INITIALIZATION_BROWN_FACTOR_HPP_
#define INITIALIZATION_BROWN_FACTOR_HPP_

#include <cmath>

#include "box.hpp"
#include "particle_parameter.hpp"
#include "time_step.hpp"

class brown_factor : virtual public box, public time_step {
 private:
  double BD_r_1_;
  double BD_r_2_;
  double BD_v_1_;
  double BD_v_2_;
  double BD_v_3_;
  double BD_g0_1_;
  double BD_g1_1_;
  double BD_g1_2_;
  double C(const double x);
  double G(const double x);

 public:
  void calc_BD_factor();
  double BD_r_1() const;
  double BD_r_2() const;
  double BD_v_1() const;
  double BD_v_2() const;
  double BD_v_3() const;
  double BD_g0_1() const;
  double BD_g1_1() const;
  double BD_g1_2() const;
  void print_brown_factor();
};

#endif  // INITIALIZATION_BROWN_FACTOR_HPP_
