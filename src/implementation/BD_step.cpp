/**
 * @file BD_step.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "BD_step.hpp"

BD_step::BD_step(/* args */) {}

void BD_step::run_BD_step() {
  std::cout << "in the BD step::run" << std::endl;
  generate_Gamma();
  calc_pos();
  calc_force();
  calc_vel();
}

BD_step::~BD_step() {}

void BD_step::calc_vel(void) {
  std::cout << "in the BD step::clac vel" << std::endl;
// calc_BD_vel();
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] =
          dr[i][ax] * BD_v_1() + f0[i][ax] * BD_v_2() + f1[i][ax] * BD_v_3();
    }
  }
}

void BD_step::calc_pos(void) {
  std::cout << "in the BD step::pos" << std::endl;
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      dr[i][ax] = v[i][ax] * BD_r_1() + f0[i][ax] * BD_r_2() + g0[i][ax];
      r[i][ax] += dr[i][ax];
      dr[i][ax] += g1[i][ax];
    }
  }
}
