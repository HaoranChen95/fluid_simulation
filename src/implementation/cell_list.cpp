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

#include "cell_list.hpp"

cell_list::cell_list(/* args */) {}

cell_list::~cell_list() {}

void cell_list::init_cell_list(const double &cell_size) {
  /**
   * @brief calculate the number and lengeth of cell
   *
   */
  for (int ax = 0; ax < 3; ax++) {
    cell_N[ax] = static_cast<uint64_t>(l_b()[ax] / cell_size);
    cell_l[ax] = l_b()[ax] / static_cast<double>(cell_N[ax]);
  }
  /**
   * @brief generate empty cell
   *
   */
  std::vector<uint64_t> empty_vector;
  std::vector<std::vector<uint64_t>> empty_vvector;
  for (int az = 0; az < cell_N[2]; az++) {
    empty_vvector.push_back(empty_vector);
  }
  std::vector<std::vector<std::vector<uint64_t>>> empty_vvvector;
  for (int ay = 0; ay < cell_N[1]; ay++) {
    empty_vvvector.push_back(empty_vvector);
  }
  for (int ax = 0; ax < cell_N[0]; ax++) {
    cell.push_back(empty_vvvector);
  }
  /**
   * @brief generate cell list of i-j paar
   *
   */
  for (uint64_t i = 0; i < Nm(); i++) {
    cell_list_ij.push_back(empty_vector);
  }
  std::cout << "initialization of cell list finished" << std::endl;
}

void cell_list::refresh_cell_list() {
  /**
   * @brief clear all the cell in the list
   *
   */
  for (uint64_t ax = 0; ax < cell_N[0]; ax++) {
    for (uint64_t ay = 0; ay < cell_N[1]; ay++) {
      for (uint64_t az = 0; az < cell_N[2]; az++) {
        cell[ax][ay][az].clear();
      }
    }
  }
  /**
   * @brief put the index in the cell list
   *
   */
  for (uint64_t i = 0; i < Nm(); i++) {
    cell[cell_of(i, 0)][cell_of(i, 1)][cell_of(i, 2)].push_back(i);
  }

  for (uint64_t i = 0; i < Nm(); i++) {
    cell_list_ij[i].clear();
    for (uint64_t ax : neighbor_i(i, 0)) {
      for (uint64_t ay : neighbor_i(i, 1)) {
        for (uint64_t az : neighbor_i(i, 2)) {
          for (uint64_t j : cell[ax][ay][az]) {
            if (j > i) {
              cell_list_ij[i].push_back(j);
            }
          }
        }
      }
    }
  }
}

uint64_t cell_list::cell_of(const uint64_t &i, const int &ax) {
  return static_cast<uint64_t>(r_in_box(i, ax) / cell_l[ax]);
}

std::array<uint64_t, 3> cell_list::neighbor_i(const uint64_t &i,
                                              const int &ax) {
  int64_t N = cell_of(i, ax);
  std::array<uint64_t, 3> resalt = {N - 1, N, N + 1};
  if (resalt[1] == 0) {
    resalt[0] = cell_N[ax] - 1;
  }
  if (resalt[2] == cell_N[ax]) {
    resalt[2] = 0;
  }

  return resalt;
}

// std::vector<uint64_t> &cell_list::neighbor_cell(const uint64_t &i) {
//   // return &(cell_list_[i]);
// }

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

// void init_parameter(void) {
//   sp.calc_BD_factor();
//   // for (int ax = 0; ax < 3; ax++) {
//   //   set_Cell_N[ax] = static_cast<uint64_t>(sp.l_b()[ax] / sp.sigma());
//   // }
// }
