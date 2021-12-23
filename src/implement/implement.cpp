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

// void cell_list(void) {
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     list[i] = -1;
//   }

//   for (int64_t cx = 0; cx < Cell_N[0] + 2; cx++) {
//     for (int64_t cy = 0; cy < Cell_N[1] + 2; cy++) {
//       for (int64_t cz = 0; cz < Cell_N[2] + 2; cz++) {
//         cell[cx][cy][cz] = -1;
//       }
//     }
//   }

//   int64_t cx, cy, cz;
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     cx = static_cast<uint64_t>(r[i][0] / Cell_l[0]);
//     cy = static_cast<uint64_t>(r[i][1] / Cell_l[1]);
//     cz = static_cast<uint64_t>(r[i][2] / Cell_l[2]);
//     list[i] = cell[cx][cy][cz];
//     cell[cx][cy][cz] = i;
//   }
// }


//   // std::vector<const ij_paar> ij_list;
//   // for (uint64_t i = 0; i < Nm; i++) {
//   //   for (uint64_t j = i + 1; j < Nm; j++) {
//   //     // ij_list.push_front({i, j});
//   //     counter++;
//   //   }
//   // }
//   // std::cout << "counter1 " << counter << std::endl;

//   // // #pragma omp parallel
//   // // {
//   // //   // #pragma omp for reduction(+ : E_pot)
//   // for (ij_paar const &ij : ij_list) {
//   //   // std::cout << "i " << ij.i << " j " << ij.j << std::endl;
//   //   E_pot += LJ(ij.i, ij.j);
//   //   counter++;
//   // }
//   // // }

//   // std::cout << "counter " << counter << std::endl;
//   // counter = 0;
// }

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
//       dr[i][ax] = v[i][ax] * sp.BD_r_1() + f0[i][ax] * sp.BD_r_2() + g0[i][ax];
//       r[i][ax] += dr[i][ax];
//       dr[i][ax] += g1[i][ax];
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
//     std::cout << "time\t" << static_cast<double>(step) * sp.h() << "\tE_kin\t"
//               << E_kin / static_cast<double>(sp.Nm()) << "\tE_pot\t"
//               << E_pot / static_cast<double>(sp.Nm()) << "\tE\t"
//               << (E_kin + E_pot) / static_cast<double>(sp.Nm()) << std::endl;
//   }
//   if (++print_E == sp.time_01()) {
//     print_E = 0;
//   }
// }

