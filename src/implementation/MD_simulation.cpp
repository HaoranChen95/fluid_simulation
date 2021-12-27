/**
 * @file MD_simulation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "MD_simulation.hpp"

MD_simulation::MD_simulation(/* args */) {}

MD_simulation::~MD_simulation() {}

void MD_simulation::MD_relaxation() {
  std::cout << "in the MD step::MD_relaxation " << half_h() << std::endl;

  for (step = 0; step < Relax_Steps(); step++) {
    calc_pos();
    calc_force();
    calc_vel();
    vel_correcter();
    print_energy();
    write_last_cfg();
  }
}

void MD_simulation::MD_implementation() {
  for (step = 0; step <= MD_Steps(); step++) {
    calc_pos();
    calc_force();
    calc_vel();
    calc_E_kin();
    print_energy();
    write_last_cfg();
    write_cfg();
    write_energy();
  }
}

void MD_simulation::calc_vel(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] += (f0[i][ax] + f1[i][ax]) * half_h();
    }
  }
}

void MD_simulation::calc_pos(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      dr[i][ax] = v[i][ax] * h() + f1[i][ax] * half_h2();
      r[i][ax] += dr[i][ax];
    }
  }
}
