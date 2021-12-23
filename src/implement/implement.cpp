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

// double minium_image(const uint64_t &i, const uint64_t &j, const int &ax) {
//   double r_ij_ax = r[j][ax] - r[i][ax];
//   r_ij_ax -=
//       sp.l_b()[ax] * floor((r_ij_ax + sp.half_l_b()[ax]) * sp.inv_l_b()[ax]);
//   return r_ij_ax;
// }

// double LJ(uint64_t i, uint64_t j) {
//   double r_ij[3], f_LJ[3];
//   double r2, r_inv2, sr_inv2, sr_inv6, sr_inv12, coeffLJ;
//   for (int ax = 0; ax < 3; ax++) {
//     r_ij[ax] = minium_image(i, j, ax);
//     r2 += r_ij[ax] * r_ij[ax];
//   }

//   if (r2 < sp.r2_cut()) {
//     r_inv2 = 1. / r2;
//     sr_inv2 = sp.sig2() * r_inv2;
//     sr_inv6 = sr_inv2 * sr_inv2 * sr_inv2;
//     sr_inv12 = sr_inv6 * sr_inv6;
//     coeffLJ = 24. * sp.epsilon() * r_inv2 * (sr_inv12 + sr_inv12 - sr_inv6);
//     for (int ax = 0; ax < 3; ax++) {
//       f_LJ[ax] = coeffLJ * r_ij[ax];
//       f1[i][ax] -= f_LJ[ax];
//       f1[j][ax] += f_LJ[ax];
//     }
//     return 4. * sp.epsilon() * (sr_inv12 - sr_inv6) + sp.epsilon();
//   }
//   return 0.;
// }

// void calc_force(void) {
// #pragma omp parallel for
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       f0[i][ax] = f1[i][ax];
//       f1[i][ax] = 0.;
//     }
//   }
//   // cell_list();

//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (uint64_t j = i + 1; j < sp.Nm(); j++) {
//       E_pot += LJ(i, j);
//     }
//   }

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

// void generate_Gamma(void) {
//   std::random_device rd{};
//   std::mt19937 gen{rd()};
//   std::normal_distribution<double> n_d(0.0, 1.);

//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       g0[i][ax] = sp.BD_g0_1() * n_d(gen);
//       g1[i][ax] = sp.BD_g1_1() * n_d(gen) + sp.BD_g1_2() * g0[i][ax];
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
