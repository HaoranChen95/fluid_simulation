/**
 * @file particle_parameter.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "particle_parameter.hpp"

void particle_parameter::m(const double input) { m_ = input; }
double particle_parameter::m() const { return m_; }

void particle_parameter::gamma(const double input) { gamma_ = input; }
double particle_parameter::gamma() const { return gamma_; }

void particle_parameter::epsilon(const double input) {
  epsilon_ = input;
}  // TODO e = kT
double particle_parameter::epsilon() const { return epsilon_; }

void particle_parameter::sigma(const double input) {
  sigma_ = input;
  sig2_ = input * input;
  r2_cut_ = input * input * pow(2., 1. / 3.);
}
double particle_parameter::sigma() const { return sigma_; }
double particle_parameter::r2_cut() const { return r2_cut_; }
double particle_parameter::sig2() const { return sig2_; }

void particle_parameter::print_particle() {
  std::cout << "====== particle parameter ======" << std::endl;
  std::cout << "m\t" << m_ << "\tgamma\t" << gamma_ << "\tsigma\t" << sigma_
            << "\tepsilon\t" << epsilon_ << "\tr2_cut\t" << r2_cut_
            << std::endl;
}