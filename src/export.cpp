/**
 * @file export.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-24
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "export.hpp"

int cfgs_pro_file = 500;
int cfgs_in_file = 0;
int cfg_file_array_Nm = 0;
std::ofstream cfg_data_file;
std::ofstream cfg_time_file;
std::ofstream energy_file;
std::ostringstream oss;

void print_cfg(void) {
  if (cfgs_in_file == 0) {
    oss << "cfg_data_pos_aV_" << cfg_file_array_Nm << "_phi_" << density
        << "_Nm_" << Nm << ".dat";
    cfg_data_file.open(oss.str(), std::ios::app);
    oss.str("");

    oss << "cfg_time_aV_" << cfg_file_array_Nm << "_phi_" << density << "_Nm_"
        << Nm << ".dat";
    cfg_time_file.open(oss.str(), std::ios::app);
    oss.str("");
  }

  if (cfgs_in_file < cfgs_pro_file) {
    for (uint64_t i = 0; i < Nm; i++) {
      cfg_data_file << r[i][0] << " " << r[i][1] << " " << r[i][2] << " "
                    << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
    }
    cfg_time_file << step * sp.h() << std::endl;
    cfgs_in_file++;
  } else {
    cfg_file_array_Nm++;
    cfgs_in_file = 0;
    cfg_data_file.close();
    cfg_time_file.close();
  }
}

void print_Energy(void) {
  // std::cout << "writing !" << std::endl;
  if (!energy_file.is_open()) {
    oss << "energy_phi_" << density << "_Nm_" << Nm << ".txt";
    energy_file.open(oss.str(), std::ios::app);
    oss.str("");
    energy_file.precision(10);
  }

  double f_sqr = 0;

  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      f_sqr += f1[i][ax] * f1[i][ax];
    }
  }

  energy_file << static_cast<double>(step) * sp.h() << " " << f_sqr << " "
              << E_kin / static_cast<double>(Nm) << " "
              << E_pot / static_cast<double>(Nm) << " "
              << (E_kin + E_pot) / static_cast<double>(Nm) << std::endl;
}
