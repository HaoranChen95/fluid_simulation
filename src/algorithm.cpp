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

#include "algorithm.hpp"

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
    cx = static_cast<uint64_t>(r1[0][i] / Cell_l[0]);
    cy = static_cast<uint64_t>(r1[1][i] / Cell_l[1]);
    cz = static_cast<uint64_t>(r1[2][i] / Cell_l[2]);
    list[i] = cell[cx][cy][cz];
    cell[cx][cy][cz] = i;
  }
}

double minium_image(const uint64_t &i, const uint64_t &j, const int &ax) {
  double r_ij_ax = r1[ax][j] - r1[ax][i];
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
      f1[ax][i] -= f_LJ[ax];
      f1[ax][j] += f_LJ[ax];
    }
    return 4. * eps * (sr_inv12 - sr_inv6) + eps;
  }
  return 0.;
}

void calc_force(void) {
  uint64_t i, j, counter{0};
  int ax, th;
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      f0[ax][i] = f1[ax][i];
      f1[ax][i] = 0;
    }
  }
  cell_list();

  // #pragma omp parallel private(i, j, th)

  // th = omp_get_thread_num();
  // #pragma omp for
  for (i = 0; i < Nm; i++) {
    for (j = i + 1; j < Nm; j++) {
      E_pot += LJ(i, j);
    }
  }
}

// void calc_fluid_force(void) {
//   calc_force();
//   uint64_t i;
//   int ax;
//   std::random_device rd{};
//   std::mt19937 gen{rd()};
//   std::normal_distribution<double> n_d(0.0, kT_2gamma_over_m);
//   for (i = 0; i < Nm; i++) {
//     for (ax = 0; ax < 3; ax++) {
//       f1[ax][i] += n_d(gen);
//     }
//   }
// }

void calc_vel(void) {
  uint64_t i;
  int ax;
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      v[ax][i] += (f0[ax][i] + f1[ax][i]) * half_dt;
    }
  }
  // for (ax = 0; ax < 3; ax++) {
  //   std::cout << " a" << ax << " " << (f0[ax][0] + f1[ax][0]) * half_dt;
  //   std::cout << " v" << ax << " " << (v[ax][0] + v[ax][0]) * half_dt;
  // }
  // std::cout << std::endl;
}

void calc_fluid_vel(void) {
  calc_vel();
  uint64_t i;
  int ax;
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      v[ax][i] =
          dr[ax][i] * const_v_1 + f0[ax][i] * const_v_2 + f1[ax][i] * const_v_3;
    }
  }
}

void calc_pos(void) {
  uint64_t i;
  int ax;
  // #ifdef _OPENMP
  // #pragma omp for
  // #endif  // _OPENMP

  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      dr[ax][i] = v[ax][i] * dt + f1[ax][i] * half_dt2;
      r1[ax][i] += dr[ax][i];
    }
  }
  // for (ax = 0; ax < 3; ax++) {
  //   std::cout<< " r " << dr[ax][0] << " v " << v[ax][0] * dt << " f "
  //             << f1[ax][0] * half_dt2;
  // }
  // std::cout << std::endl;
}

void calc_fluid_pos(void) {
  uint64_t i;
  int ax;
  // #ifdef _OPENMP
  // #pragma omp for
  // #endif  // _OPENMP

  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      dr[ax][i] = v[ax][i] * const_r_1 + f0[ax][i] * const_r_2 + g0[ax][i];
      r1[ax][i] += dr[ax][i];
      dr[ax][i] += g1[ax][i];
    }
  }
}

void calc_E_kin(void) {
  E_kin = 0.;
  uint64_t i;
  int ax;
  // th = omp_get_thread_num();
  // #pragma omp for reduction(+ : E_kin)
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      E_kin += v[ax][i] * v[ax][i];
    }
  }
  E_kin *= 0.5;
}

double calc_f_all(void) {
  double f_all = 0;
  uint64_t i;
  int ax;
  // #pragma omp for reduction(+ : E_kin)
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      f_all += v[ax][i];
    }
  }
  return f_all;
}

void MD_Step(void) {
  E_pot = 0;
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
  std::cout << " E_kin " << E_kin / static_cast<double>(Nm) << " E_pot "
            << E_pot / static_cast<double>(Nm) << " E "
            << (E_kin + E_pot) / static_cast<double>(Nm) << " half dt "
            << half_dt
            // << " v_all " << calc_f_all()
            << std::endl;
}
