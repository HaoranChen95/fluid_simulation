/**
 * @file algorithm.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-18
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "implement.hpp"

void cell_list(void) {
  for (uint64_t i = 0; i < Nm; i++) {
    list[i] = -1;
  }

  for (int64_t cx = 0; cx < Cell_N[0] + 2; cx++) {
    for (int64_t cy = 0; cy < Cell_N[1] + 2; cy++) {
      for (int64_t cz = 0; cz < Cell_N[2] + 2; cz++) {
        cell[cx][cy][cz] = -1;
      }
    }
  }

  int64_t cx, cy, cz;
  for (uint64_t i = 0; i < Nm; i++) {
    cx = static_cast<uint64_t>(r[i][0] / Cell_l[0]);
    cy = static_cast<uint64_t>(r[i][1] / Cell_l[1]);
    cz = static_cast<uint64_t>(r[i][2] / Cell_l[2]);
    list[i] = cell[cx][cy][cz];
    cell[cx][cy][cz] = i;
  }
}

double minium_image(const uint64_t &i, const uint64_t &j, const int &ax) {
  double r_ij_ax = r[j][ax] - r[i][ax];
  r_ij_ax -= l_b[ax] * floor((r_ij_ax + half_l_b[ax]) * inv_l_b[ax]);
  return r_ij_ax;
}

double LJ(uint64_t i, uint64_t j) {
  double r_ij[3], f_LJ[3];
  double r2, r_inv2, sr_inv2, sr_inv6, sr_inv12, coeffLJ;
  for (int ax = 0; ax < 3; ax++) {
    r_ij[ax] = minium_image(i, j, ax);
    r2 += r_ij[ax] * r_ij[ax];
  }

  if (r2 < r2_cut) {
    r_inv2 = 1. / r2;
    sr_inv2 = sig2 * r_inv2;
    sr_inv6 = sr_inv2 * sr_inv2 * sr_inv2;
    sr_inv12 = sr_inv6 * sr_inv6;
    coeffLJ = 24. * eps * r_inv2 * (sr_inv12 + sr_inv12 - sr_inv6);
    for (int ax = 0; ax < 3; ax++) {
      f_LJ[ax] = coeffLJ * r_ij[ax];
      f1[i][ax] -= f_LJ[ax];
      f1[j][ax] += f_LJ[ax];
    }
    return 4. * eps * (sr_inv12 - sr_inv6) + eps;
  }
  return 0.;
}

void calc_force(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      f0[i][ax] = f1[i][ax];
      f1[i][ax] = 0.;
    }
  }
  // cell_list();

  // uint64_t counter{0};

  for (uint64_t i = 0; i < Nm; i++) {
    for (uint64_t j = i + 1; j < Nm; j++) {
      E_pot += LJ(i, j);
      // counter++;
    }
  }

  // std::vector<const ij_paar> ij_list;
  // for (uint64_t i = 0; i < Nm; i++) {
  //   for (uint64_t j = i + 1; j < Nm; j++) {
  //     // ij_list.push_front({i, j});
  //     counter++;
  //   }
  // }
  // std::cout << "counter1 " << counter << std::endl;

  // // #pragma omp parallel
  // // {
  // //   // #pragma omp for reduction(+ : E_pot)
  // for (ij_paar const &ij : ij_list) {
  //   // std::cout << "i " << ij.i << " j " << ij.j << std::endl;
  //   E_pot += LJ(ij.i, ij.j);
  //   counter++;
  // }
  // // }

  // std::cout << "counter " << counter << std::endl;
  // counter = 0;
}

void calc_vel(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] += (f0[i][ax] + f1[i][ax]) * sp.half_h();
    }
  }
}

void calc_fluid_vel(void) {
  calc_vel();
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] =
          dr[i][ax] * const_v_1 + f0[i][ax] * const_v_2 + f1[i][ax] * const_v_3;
    }
  }
}

void calc_pos(void) {
#pragma omp parallel for
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      dr[i][ax] = v[i][ax] * sp.h() + f1[i][ax] * sp.half_h2();
      r[i][ax] += dr[i][ax];
    }
  }
}

void calc_fluid_pos(void) {
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      dr[i][ax] = v[i][ax] * const_r_1 + f0[i][ax] * const_r_2 + g0[i][ax];
      r[i][ax] += dr[i][ax];
      dr[i][ax] += g1[i][ax];
    }
  }
}

void calc_E_kin(void) {
  E_kin = 0.;
#pragma omp for reduction(+ : E_kin)
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      E_kin += v[i][ax] * v[i][ax];
    }
  }
  E_kin *= 0.5;
}

void MD_Step(void) {
  E_pot = 0.;
  if (open_fluid) {
    init_gamma();
    calc_fluid_pos();
    calc_force();
    calc_fluid_vel();
  } else {

    calc_pos();
    calc_force();
    calc_vel();
  }

  calc_E_kin();
  if (print_E == 0) {
    std::cout << "time\t" << static_cast<double>(step) * sp.h() << "\tE_kin\t"
              << E_kin / static_cast<double>(Nm) << "\tE_pot\t"
              << E_pot / static_cast<double>(Nm) << "\tE\t"
              << (E_kin + E_pot) / static_cast<double>(Nm)
              // << " v_all " << calc_f_all()
              << std::endl;
  }
  if (++print_E == sp.time_01()) {
    print_E = 0;
  }
}

void vel_correcter(void) {
  calc_E_kin();
  double a = sqrt(1.5 * sp.kT() * Nm / E_kin);
  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      v[i][ax] *= a;
    }
  }
}
