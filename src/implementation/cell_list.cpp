/**
 * @file cell_list.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief 
 * @version 0.1
 * @date 2021-12-26
 * 
 * @copyright Copyright (c) 2021
 * 
 */

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
// }


// uint64_t set_Cell_N[3];
// const uint64_t *const Cell_N = set_Cell_N;

// double set_Cell_l[3];
// const double *const Cell_l = set_Cell_l;

// int64_t ***cell;
// int64_t *list;

// double E_kin, E_pot;
// int print_E;
// int FREQ_print_E = 10;

// void init_parameter(void) {
//   sp.calc_BD_factor();
//   // for (int ax = 0; ax < 3; ax++) {
//   //   set_Cell_N[ax] = static_cast<uint64_t>(sp.l_b()[ax] / sp.sigma());
//   // }
// }