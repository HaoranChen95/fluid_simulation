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
  fluid_simulation fluid(argc, argv);
  fluid.relax();
  fluid.implement();
  // fluid.implementation();

  // int FREQ_CFG_DETA = 1;
  // init_system(argc, argv);
  // calc_force();
  // std::cout << "there !!! " << sp.BD_r_1() << " " << sp.BD_r_2() << " "
  //           << sp.BD_v_1() << " " << sp.BD_v_2() << " " << sp.BD_v_3() << " "
  //           << sp.BD_g0_1() << " " << sp.BD_g1_1() << " " << sp.BD_g1_2() << std::endl;

  // sp.Relax_Steps(0);
  // for (step = 0; step < sp.Relax_Steps(); step++) {
  //   MD_Step();
  //   if (!sp.gamma()) {
  //     vel_correcter();
  //   }
  // }
  // for (step = 0; step < sp.MD_Steps(); step++) {
  //   MD_Step();
  //   if (step % sp.time_01() == 0) {
  //     print_Energy();
  //   }
  //   if (step % sp.time_01() == 0) {
  //     print_cfg();
  //   }
  // }

  // close_system();

  return 0;
}
