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

// void calc_BD_pos(void) {
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       dr[i][ax] = v[i][ax] * sp.BD_r_1() + f0[i][ax] * sp.BD_r_2() +
//       g0[i][ax]; r[i][ax] += dr[i][ax]; dr[i][ax] += g1[i][ax];
//     }
//   }
// }
