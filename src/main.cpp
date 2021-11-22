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
  // Relax_Steps = 100000;
  // MD_Steps = 10000000;
  Relax_Steps = 0;
  MD_time = 1.;
  set_dt = 0.001;
  MD_Steps = static_cast<uint64_t>(MD_time / dt);

  int FREQ_CFG_DETA = 1;
  read_arg(argc, argv);
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
    // std::cout << "there !!! " << step << std::endl;

    if (step % time_01 == 0) {
      print_Energy();
    }
    if (step % time_01 == 0) {
      print_cfg();
    }
  }

  close_system();

  return 0;
}
