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

#include "output.hpp"

output::output(/* args */) {}

output::~output() {}

uint64_t output::cfg_freq() {
  if (step < time_1000()) {
    int step_remain = step % time_10();
    if (step_remain == 0) {
      cfg_freq_ = time_0001();
    }
    if (step_remain == time_01()) {
      cfg_freq_ = time_001();
    }
    if (step_remain == time_1()) {
      cfg_freq_ = time_01();
    }
    // if (step_remain == time_10()) {
    //   cfg_freq_ = time_1();
    // }
  } else if (step == time_1000()) {
    cfg_freq_ = time_1();
  }
  return cfg_freq_;
}

void output::print_energy() {
  if (step % time_01() == 0) {
    std::cout << "time\t" << step_time() << "\tE_pot\t" << E_pot()
              << "\tE_kin\t" << E_kin() << "\tE_sum\t" << E_kin() + E_pot()
              << std::endl;
  }
}

void output::write_last_cfg(void) {
  if (step % time_1() == 0) {
    oss << "read_init_cfg_vel_Nm_" << Nm() << "_kT_" << kT() << ".txt";
    fn = oss.str();
    last_cfg_pos_vel.open(fn, std::ios::trunc);
    oss.str("");
    last_cfg_pos_vel.precision(6);
    for (int i = 0; i < Nm(); ++i) {
      last_cfg_pos_vel << r[i][0] << " " << r[i][1] << " " << r[i][2] << " "
                       << v[i][0] << " " << v[i][1] << " " << v[i][2]
                       << std::endl;
    }
    last_cfg_pos_vel.close();
  }
}

void output::write_cfg(void) {
  if (step % cfg_freq() == 0) {
    if (cfgs_in_file == 0) {
      oss << "cfg_data_aV_" << cfg_file_array_Nm << "_phi_" << density()
          << "_Nm_" << Nm() << "_kT_" << kT() << "_gam_" << gamma() << ".txt";
      cfg_data_file.open(oss.str(), std::ios::app);
      oss.str("");

      oss << "cfg_time_aV_" << cfg_file_array_Nm << "_phi_" << density()
          << "_Nm_" << Nm() << "_kT_" << kT() << "_gam_" << gamma() << ".txt";
      cfg_time_file.open(oss.str(), std::ios::app);
      oss.str("");
    }

    if (cfgs_in_file < cfgs_pro_file) {
      for (uint64_t i = 0; i < Nm(); i++) {
        cfg_data_file << r[i][0] << " " << r[i][1] << " " << r[i][2] << " "
                      << v[i][0] << " " << v[i][1] << " " << v[i][2]
                      << std::endl;
      }
      cfg_time_file << step * h() << std::endl;
      cfgs_in_file++;
    } else {
      cfg_file_array_Nm++;
      cfgs_in_file = 0;
      cfg_data_file.close();
      cfg_time_file.close();
    }
  }
}

void output::write_energy(void) {
  if (step % time_01() == 0) {
    if (!energy_file.is_open()) {
      oss << "energy_phi_" << density() << "_Nm_" << Nm() << "_kT_" << kT()
          << "_gam_" << gamma() << ".txt";
      energy_file.open(oss.str(), std::ios::app);
      oss.str("");
      energy_file.precision(10);
    }

    energy_file << step_time() << " " << E_kin() << " " << E_pot() << " "
                << E_kin() + E_pot() << std::endl;
  }
}
