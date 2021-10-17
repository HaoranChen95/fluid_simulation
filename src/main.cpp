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
  // Relax_Steps = 100000;
  // MD_Steps = 10000000;
  Relax_Steps = 0;
  MD_Steps = 100000;
  set_dt = 0.001;

  int FREQ_CFG_DETA = 1;

  read_config();
  init_system();
  calc_force();

  std::cout << "half lx = " << half_l_b[0] << std::endl
            << "gamma = " << gam << std::endl;
  for (step = 0; step < Relax_Steps; step++) {
    MD_Step();
    if (!open_fluid) {
      vel_correcter();
    }
  }
  for (step = 0; step < MD_Steps; step++) {
    MD_Step();

    if (step == 100) {
      FREQ_CFG_DETA = 10;
    }
    if (step == 1000) {
      FREQ_CFG_DETA = 100;
    }
    if (step == 10000) {
      FREQ_CFG_DETA = 1000;
    }
    if (step == 100000) {
      FREQ_CFG_DETA = 10000;
    }
    if (step == 1000000) {
      FREQ_CFG_DETA = 100000;
    }

    if (step % FREQ_CFG_DETA == 0) {
      print_cfg();
    }
  }

  close_system();

  return 0;
}
