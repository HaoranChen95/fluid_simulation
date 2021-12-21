/**
 * @file particle_parameter.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef PARTICLE_PARAMETER_HPP_
#define PARTICLE_PARAMETER_HPP_

#include <cmath>
#include <iostream>

class particle_parameter {
 private:
  double m_;
  double gamma_;
  double sigma_;
  double epsilon_;
  double r2_cut_;
  double sig2_;

 public:
  void m(const double input);
  double m() const;

  void gamma(const double input);
  double gamma() const;
  void epsilon(const double input);
  double epsilon() const;

  void sigma(const double input);
  double sigma() const;
  double r2_cut() const;
  double sig2() const;
};

#endif  // PARTICLE_PARAMETER_HPP_