/**
 * @file velocity.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "velocity.hpp"

void velocity::init_velocity() {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, kT());

  std::array<double, 3> v_sum = {0, 0, 0};

  std::array<double, 3> new_v;
  for (uint64_t i = 0; i < Nm(); i++) {
    new_v[0] = n_d(gen);
    new_v[1] = n_d(gen);
    new_v[2] = n_d(gen);
    v.push_back(new_v);
  }

  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      v_sum[ax] += v[i][ax];
    }
  }

  for (int ax = 0; ax < 3; ax++) {
    v_sum[ax] /= static_cast<double>(Nm());
  }

  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] -= v_sum[ax];
    }
  }

  std::cout << "initialization of velocity finished: size " << v.size()
            << std::endl;
}

void velocity::calc_E_kin() {
  E_kin_ = 0.;
#pragma omp for reduction(+ : E_kin_)
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      E_kin_ += v[i][ax] * v[i][ax];
    }
  }
  E_kin_ *= 0.5 / static_cast<double>(Nm());
}

void velocity::calc_mean_E_kin() {
  step_++;
  mean_E_kin_ = (mean_E_kin_ * (static_cast<double>(step_ - 1)) + E_kin_) /
                static_cast<double>(step_);
}

double velocity::E_kin() const { return E_kin_; }

void velocity::vel_correcter(void) {
  calc_E_kin();
  calc_mean_E_kin();
  double a = sqrt(1.5 * kT() / mean_E_kin_);
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] *= a;
    }
  }
}

velocity::~velocity() {}
