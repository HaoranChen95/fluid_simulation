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

#ifndef INITIALIZATION_PARTICLE_PARAMETER_HPP_
#define INITIALIZATION_PARTICLE_PARAMETER_HPP_

#ifdef _OPENMP
#define N_THREADS 1

#else
#define N_THREADS 1
#endif  // _OPENMP

#include <omp.h>

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
  double sig3_;

 public:
  void m(const double input);
  void gamma(const double input);
  void epsilon(const double input);
  void sigma(const double input);

  double m() const;
  double gamma() const;
  double epsilon() const;
  double sigma() const;

  double r2_cut() const;
  double sig2() const;
  double sig3() const;
  void print_particle();
};

#endif  // INITIALIZATION_PARTICLE_PARAMETER_HPP_
