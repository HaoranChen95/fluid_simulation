/**
 * @file BD_simulation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "BD_simulation.hpp"

BD_simulation::BD_simulation() {}

void BD_simulation::BD_relaxation() {
  std::cout << "in the BD step::MD_relaxation " << std::endl;

  for (step = 0; step < Relax_Steps(); step++) {
    generate_Gamma();
    calc_pos();
    calc_force();
    calc_vel();
    calc_E_kin();
    print_energy();
    write_last_cfg();
  }
}

void BD_simulation::BD_implementation() {
  for (step = 0; step <= MD_Steps(); step++) {
    generate_Gamma();
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

BD_simulation::~BD_simulation() {}

void BD_simulation::calc_vel(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] =
          dr[i][ax] * BD_v_1() + f0[i][ax] * BD_v_2() + f1[i][ax] * BD_v_3();
    }
  }
}

void BD_simulation::calc_pos(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      dr[i][ax] = v[i][ax] * BD_r_1() + f0[i][ax] * BD_r_2() + g0[i][ax];
      r[i][ax] += dr[i][ax];
      dr[i][ax] += g1[i][ax];
    }
  }
}
