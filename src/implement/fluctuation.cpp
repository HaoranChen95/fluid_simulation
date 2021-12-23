/**
 * @file fluctuation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief 
 * @version 0.1
 * @date 2021-12-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "fluctuation.hpp"


// void init_gamma(void) {
//   std::random_device rd{};
//   std::mt19937 gen{rd()};
//   std::normal_distribution<double> n_d(0.0, 1.);

//   std::array<double, 3> new_g0, new_g1;

//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       new_g0[ax] = sp.BD_g0_1() * n_d(gen);
//       new_g1[ax] = sp.BD_g1_1() * n_d(gen) + sp.BD_g1_2() * new_g0[ax];
//     }
//     g0.push_back(new_g0);
//     g1.push_back(new_g1);
//   }
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