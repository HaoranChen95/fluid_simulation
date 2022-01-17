/**
 * @file cell_list.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef INITIALIZATION_CELL_LIST_HPP_
#define INITIALIZATION_CELL_LIST_HPP_

#include <array>
#include <vector>

#include "position.hpp"

class cell_list : virtual protected position {
 private:
  std::vector<std::vector<std::vector<std::vector<uint64_t>>>> cell;
  std::array<uint64_t, 3> cell_N;
  std::array<double, 3> cell_l;
  uint64_t cell_of(const uint64_t &i, const int &ax);
  std::array<uint64_t, 3> neighbor_i(const uint64_t &i, const int &ax);

 protected:
  std::vector<std::vector<uint64_t>> cell_list_ij;

 public:
  void init_cell_list(const double &cell_size);
  void refresh_cell_list();
  //   std::vector<uint64_t> &neighbor_cell(const uint64_t &i);
};

#endif  // INITIALIZATION_CELL_LIST_HPP_
