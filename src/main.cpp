/**
 * @file main.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-16
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "main.hpp"

int main(const int argc, const char* argv[]) {
  std::cout << "Hello" << std::endl;

  MD_Steps = 1000;
  set_dt = 0.0001;

  read_config();
  init_system();
  calc_force();

  std::cout << "half lx = " << half_l_b[0] << std::endl
            << "gamma = " << gam << std::endl;
  for (step = 0; step < Relax_Steps; step++) {
    MD_Step();
  }
  for (step = 0; step < MD_Steps; step++) {
    MD_Step();
    print_cfg();
  }

  std::cout << "kT 2 gamma over m " << kT_2gamma_over_m << std::endl;
  // r1[0][0] = 0;
  // r1[0][1] = 151;
  std::cout << "const 1 = " << const_r_1 << std::endl;
  std::cout << "const 2 = " << const_r_2 << std::endl;
  std::cout << "const 3 = " << const_v_1 << std::endl;
  std::cout << "const 4 = " << const_v_2 << std::endl;
  std::cout << "const 5 = " << const_v_3 << std::endl;
  std::cout << "const 6 = " << const_g0_1 << std::endl;
  std::cout << "const 7 = " << const_g1_1 << std::endl;
  std::cout << "const 8 = " << const_g1_2 << std::endl;
  // r1[0][1] = 101;
  // std::cout << "r ij 2 = " << minium_image(0, 1, 0) << std::endl;
  // r1[0][1] = 51;
  // std::cout << "r ij 3 = " << minium_image(0, 1, 0) << std::endl;
  // r1[0][1] = 1;
  // std::cout << "r ij 4 = " << minium_image(0, 1, 0) << std::endl;
  // r1[0][1] = -49;
  // std::cout << "r ij 5 = " << minium_image(0, 1, 0) << std::endl;
  // r1[0][1] = -99;
  // std::cout << "r ij 6 = " << minium_image(0, 1, 0) << std::endl;
  // r1[0][1] = -149;
  // std::cout << "r ij 7 = " << minium_image(0, 1, 0) << std::endl;

  close_system();

  return 0;
}
