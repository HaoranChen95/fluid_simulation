/**
 * @file position.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "position.hpp"

// double position::r_in_box(const uint64_t &i, const int &ax) {
//   return r[i][ax] - l_b()[ax] * floor(r[i][ax] * inv_l_b()[ax]);
// }

void position::init_position() {
  int row_x, row_y, row_z;
  double root2 = sqrt(2.);
  double half_root2 = sqrt(0.5);

  row_x = static_cast<int>(l_b()[0] / root2 / sigma());
  row_y = static_cast<int>(l_b()[1] / root2 / sigma());
  row_z = static_cast<int>(l_b()[2] / root2 / sigma());

  uint64_t i = 0;
  std::array<double, 3> new_r;
  for (int r_z = 0; r_z < row_z; r_z++) {
    for (int r_y = 0; r_y < row_y; r_y++) {
      for (int r_x = 0; r_x < row_x; r_x++) {
        new_r[0] = r_x * sigma() * root2;
        new_r[1] = r_y * sigma() * root2;
        new_r[2] = r_z * sigma() * root2;
        r.push_back(new_r);
        dr.push_back({0, 0, 0});
        if (++i >= Nm()) {
          goto finish;
        }

        new_r[0] = r_x * sigma() * root2 + half_root2;
        new_r[1] = r_y * sigma() * root2 + half_root2;
        new_r[2] = r_z * sigma() * root2;
        r.push_back(new_r);
        dr.push_back({0, 0, 0});
        if (++i >= Nm()) {
          goto finish;
        }

        new_r[0] = r_x * sigma() * root2 + half_root2;
        new_r[1] = r_y * sigma() * root2;
        new_r[2] = r_z * sigma() * root2 + half_root2;
        r.push_back(new_r);
        dr.push_back({0, 0, 0});
        if (++i >= Nm()) {
          goto finish;
        }

        new_r[0] = r_x * sigma() * root2;
        new_r[1] = r_y * sigma() * root2 + half_root2;
        new_r[2] = r_z * sigma() * root2 + half_root2;
        r.push_back(new_r);
        dr.push_back({0, 0, 0});
        if (++i >= Nm()) {
          goto finish;
        }
      }
    }
  }
  Nm(r.size());
  calc_density();
finish:
  std::cout << "initialization of position finished: size " << r.size()
            << std::endl;
}

double position::minium_image(const uint64_t &i, const uint64_t &j,
                              const int &ax) {
  double r_ij_ax = r[j][ax] - r[i][ax];
  r_ij_ax -= l_b()[ax] * floor((r_ij_ax + half_l_b()[ax]) * inv_l_b()[ax]);
  return r_ij_ax;
}

position::~position() {}
