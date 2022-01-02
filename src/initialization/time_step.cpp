/**
 * @file time_step.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-05
 *
 * @copyright Copyright (c) 2021
 *
 */
#include "time_step.hpp"

void time_step::Relax_time(const uint64_t input) {
  Relax_time_ = input;
  if (h_ > 0) {
    Relax_Steps_ = static_cast<uint64_t>(input / h_);
  }
}

void time_step::MD_time(const double input) {
  MD_time_ = input;
  if (h_ > 0) {
    MD_Steps_ = static_cast<uint64_t>(MD_time_ / h_);
  }
}

void time_step::h(const double input) {
  h_ = input;
  half_h_ = 0.5 * h_;
  half_h2_ = 0.5 * h_ * h_;
  time_1000_ = static_cast<uint64_t>(1000. / h_);
  time_100_ = static_cast<uint64_t>(100. / h_);
  time_10_ = static_cast<uint64_t>(10. / h_);
  time_1_ = static_cast<uint64_t>(1. / h_);
  time_01_ = static_cast<uint64_t>(0.1 / h_);
  time_001_ = static_cast<uint64_t>(0.01 / h_);
  time_0001_ = static_cast<uint64_t>(0.001 / h_);
  MD_Steps_ = static_cast<uint64_t>(MD_time_ / h_);
}
uint64_t time_step::MD_Steps() const { return MD_Steps_; }
uint64_t time_step::Relax_Steps() const { return Relax_Steps_; }
double time_step::step_time() const { return static_cast<double>(step) * h_; }

/**
 * @brief out put time step h
 *
 * @return double
 */
double time_step::h() const { return h_; }
/**
 * @brief out put half time step
 *
 * @return double
 */
double time_step::half_h() const { return half_h_; }
double time_step::half_h2() const { return half_h2_; }
uint64_t time_step::time_0001() const { return time_0001_; }
uint64_t time_step::time_001() const { return time_001_; }
uint64_t time_step::time_01() const { return time_01_; }
uint64_t time_step::time_1() const { return time_1_; }
uint64_t time_step::time_10() const { return time_10_; }
uint64_t time_step::time_100() const { return time_100_; }
uint64_t time_step::time_1000() const { return time_1000_; }

void time_step::print_time() {
  std::cout << "====== time parameter ======" << std::endl;
  std::cout << "time length\t" << h_ << "\trelax time\t" << Relax_time_
            << "\trun time\t" << MD_time_ << std::endl;
}
