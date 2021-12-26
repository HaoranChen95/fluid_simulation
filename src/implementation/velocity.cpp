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

void velocity::init_velocity(const uint64_t Nm, const double kT) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, kT);
  uint64_t i;
  int ax, th;
  double v_sum[3] = {0, 0, 0};

  std::array<double, 3> new_v;
  for (i = 0; i < Nm; i++) {
    new_v[0] = n_d(gen);
    new_v[1] = n_d(gen);
    new_v[2] = n_d(gen);
    v.push_back(new_v);
  }
#pragma omp parallel for
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      v_sum[ax] += v[i][ax];
    }
  }
  std::cout << "initialization of velocity finished" << std::endl;
}

// void calc_E_kin(void) {
//   E_kin = 0.;
// // double v2_max = 0;
// // uint64_t i_max = 0;
// // double v2 = 0;
// #pragma omp for reduction(+ : E_kin)
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     // v2 = 0;
//     for (int ax = 0; ax < 3; ax++) {
//       E_kin += v[i][ax] * v[i][ax];
//     }
//     // if (v2 > v2_max) {
//     //   v2_max = v2;
//     //   i_max = i;
//     // }
//   }
//   // std::cout << "i_max " << i_max << " v2_max " << v2_max << " dr "
//   //           << dr[i_max][0] << " " << dr[i_max][1] << " " << dr[i_max][2]
//   //           << " vdt " << v[i_max][0] * sp.h() << " " << v[i_max][1] * sp.h()
//   //           << " " << v[i_max][2] * sp.h() << std::endl
//   //           << " g " << g0[i_max][0] << " " << g0[i_max][1] << " "
//   //           << g0[i_max][2] << " " << g1[i_max][0] << " " << g1[i_max][1] <<
//   //           " "
//   //           << g1[i_max][2] << std::endl;
//   E_kin *= 0.5;
// }


// void vel_correcter(void) {
//   calc_E_kin();
//   double a = sqrt(1.5 * sp.kT() * sp.Nm() / E_kin);
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       v[i][ax] *= a;
//     }
//   }
// }

velocity::~velocity() {}
