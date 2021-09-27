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
std::ostringstream oss;

void print_cfg(void) {
  if (cfgs_in_file == 0) {
    oss << "cfg_data_pos_aV_" << cfg_file_array_Nm << "_Nm_" << Nm  << ".txt";
    cfg_data_file.open(oss.str(), std::ios::app);
    oss.str("");

    oss << "cfg_time_aV_" << cfg_file_array_Nm << "_Nm_" << Nm  << ".txt";
    cfg_time_file.open(oss.str(), std::ios::app);
    oss.str("");
  }

  if (cfgs_in_file < cfgs_pro_file) {
    for (uint64_t i = 0; i < Nm; i++) {
      cfg_data_file << r1[0][i] << " " << r1[1][i] << " " << r1[2][i]
                    << std::endl;
    }
    cfg_time_file << step * dt << std::endl;
    cfgs_in_file++;
  } else {
    std::cout << "set new file" << std::endl 
    << std::endl << std::endl 
    << std::endl << std::endl;
    cfg_file_array_Nm++;
    cfgs_in_file = 0;
    cfg_data_file.close();
    cfg_time_file.close();
  }
}
