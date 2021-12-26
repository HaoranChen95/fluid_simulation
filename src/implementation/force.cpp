/**
 * @file force.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "force.hpp"

void force::init_force() {
  std::array<double, 3> new_f{0, 0, 0};
  for (uint64_t i = 0; i < Nm(); i++) {
    f0.push_back(new_f);
    f1.push_back(new_f);
  }
  calc_force();
  std::cout << "initialization of force finished" << std::endl;
}

void force::calc_force() {
  E_pot_ = 0;
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      f0[i][ax] = f1[i][ax];
      f1[i][ax] = 0.;
    }
  }
  // cell_list();

#pragma omp for reduction(+ : E_pot_)
  for (uint64_t i = 0; i < Nm(); i++) {
    for (uint64_t j = i + 1; j < Nm(); j++) {
      E_pot_ += LJ(i, j);
    }
  }
  E_pot_ /= static_cast<double>(Nm());
}

double force::E_pot() const { return E_pot_; }

double force::LJ(uint64_t i, uint64_t j) {
  double r_ij[3], f_LJ[3];
  double r2, r_inv2, sr_inv2, sr_inv6, sr_inv12, coeffLJ;
  for (int ax = 0; ax < 3; ax++) {
    r_ij[ax] = minium_image(i, j, ax);
    r2 += r_ij[ax] * r_ij[ax];
  }

  if (r2 < r2_cut()) {
    r_inv2 = 1. / r2;
    sr_inv2 = sig2() * r_inv2;
    sr_inv6 = sr_inv2 * sr_inv2 * sr_inv2;
    sr_inv12 = sr_inv6 * sr_inv6;
    coeffLJ = 24. * epsilon() * r_inv2 * (sr_inv12 + sr_inv12 - sr_inv6);
    for (int ax = 0; ax < 3; ax++) {
      f_LJ[ax] = coeffLJ * r_ij[ax];
      f1[i][ax] -= f_LJ[ax];
      f1[j][ax] += f_LJ[ax];
    }
    return 4. * epsilon() * (sr_inv12 - sr_inv6) + epsilon();
  }
  return 0.;
}

force::~force() {}
