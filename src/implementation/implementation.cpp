/**
 * @file implementation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "implementation.hpp"

implementation::implementation(initialization system) {
  std::cout << "in the implementation" << std::endl;
}

// void calc_MD_vel(void) {
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       v[i][ax] += (f0[i][ax] + f1[i][ax]) * sp.half_h();
//     }
//   }
// }

// void calc_BD_vel(void) {
// // calc_MD_vel();
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       v[i][ax] = dr[i][ax] * sp.BD_v_1() + f0[i][ax] * sp.BD_v_2() +
//                  f1[i][ax] * sp.BD_v_3();
//     }
//   }
// }

// void calc_MD_pos(void) {
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       dr[i][ax] = v[i][ax] * sp.h() + f1[i][ax] * sp.half_h2();
//       r[i][ax] += dr[i][ax];
//     }
//   }
// }

// void calc_BD_pos(void) {
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       dr[i][ax] = v[i][ax] * sp.BD_r_1() + f0[i][ax] * sp.BD_r_2() +
//       g0[i][ax]; r[i][ax] += dr[i][ax]; dr[i][ax] += g1[i][ax];
//     }
//   }
// }

// void MD_Step(void) {
//   E_pot = 0.;
//   if (sp.gamma()) {
//     generate_Gamma();
//     calc_BD_pos();
//     // calc_MD_pos();
//     calc_force();
//     calc_BD_vel();
//   } else {
//     calc_MD_pos();
//     calc_force();
//     calc_MD_vel();
//   }

//   calc_E_kin();
//   if (print_E == 0) {
//     std::cout << "time\t" << static_cast<double>(step) * sp.h() <<
//     "\tE_kin\t"
//               << E_kin / static_cast<double>(sp.Nm()) << "\tE_pot\t"
//               << E_pot / static_cast<double>(sp.Nm()) << "\tE\t"
//               << (E_kin + E_pot) / static_cast<double>(sp.Nm()) << std::endl;
//   }
//   if (++print_E == sp.time_01()) {
//     print_E = 0;
//   }
// }

implementation::~implementation() {}
