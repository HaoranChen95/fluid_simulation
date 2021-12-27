/**
 * @file time_step.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-05
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef INITIALIZATION_TIME_STEP_HPP_
#define INITIALIZATION_TIME_STEP_HPP_

#include <cstdint>

class time_step {
 private:
  double MD_time_;
  double Relax_time_;
  double h_;
  double half_h_;
  double half_h2_;
  uint64_t MD_Steps_;
  uint64_t Relax_Steps_;
  uint64_t time_0001_;
  uint64_t time_001_;
  uint64_t time_01_;
  uint64_t time_1_;
  uint64_t time_10_;
  uint64_t time_100_;
  uint64_t time_1000_;

 public:
  uint64_t step;

  void MD_Steps(const double input);
  void Relax_time(const uint64_t input);
  void MD_time(const double input);
  void h(const double input);

  uint64_t MD_Steps() const;
  uint64_t Relax_Steps() const;
  double MD_time() const;
  double step_time() const;
  double h() const;
  double half_h() const;
  double half_h2() const;
  uint64_t time_0001() const;
  uint64_t time_001() const;
  uint64_t time_01() const;
  uint64_t time_1() const;
  uint64_t time_10() const;
  uint64_t time_100() const;
  uint64_t time_1000() const;
};

#endif  // INITIALIZATION_TIME_STEP_HPP_
