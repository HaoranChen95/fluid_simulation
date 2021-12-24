/**
 * @file init.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-16
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "initialization.hpp"

initialization::initialization(const int argc, const char **argv) {
  read_arg(argc, argv);
  read_config();
  init_position();
  init_velocity(Nm(), kT());
  init_force();
  if (gamma()) {
    init_fluctuation();
  }
  std::cout << "v " << v.size() << std::endl;
  std::cout << "r " << r.size() << std::endl;
  std::cout << "f " << f0.size() << std::endl;
  std::cout << "E_pot " << E_pot() << std::endl;
}

void initialization::read_arg(const int argc, const char **argv) {
  MD_time(std::stod(argv[1]));
  h(std::stod(argv[2]));
  density(std::stod(argv[3]));
  gamma(std::stod(argv[4]));
}

void initialization::read_config() {
  std::cout << "read initial file:" << std::endl;

  std::ifstream config_file;
  std::string line;
  // std::string head;
  int split_at;
  config_file.open("config.txt", std::ios::in);
  if (config_file.is_open()) {
    while (std::getline(config_file, line)) {
      if (line[0] != '#') {
        split_at = line.find('=');
        std::cout << line.substr(0, split_at - 1) << ":\t"
                  << line.substr(split_at + 2) << std::endl;

        const std::string &head = line.substr(0, split_at - 1);
        const std::string &value = line.substr(split_at + 2);

        if (head == "kT") {
          kT(std::stod(value));
        } else if (head == "lx") {
          l_b(0, std::stod(value));
        } else if (head == "ly") {
          l_b(1, std::stod(value));
        } else if (head == "lz") {
          l_b(2, std::stod(value));
        } else if (head == "sigma") {
          sigma(std::stod(value));
        } else if (head == "epsilon") {
          epsilon(std::stod(value));
        } else if (head == "m") {
          m(std::stod(value));
        }
      }
    }
    calc_Nm();
    calc_density();
  }
}

// int64_t step;

// uint64_t set_Cell_N[3];
// const uint64_t *const Cell_N = set_Cell_N;

// double set_Cell_l[3];
// const double *const Cell_l = set_Cell_l;

// int64_t ***cell;
// int64_t *list;

// double E_kin, E_pot;
// int print_E;
// int FREQ_print_E = 10;

// void init_parameter(void) {
//   sp.calc_BD_factor();
//   // for (int ax = 0; ax < 3; ax++) {
//   //   set_Cell_N[ax] = static_cast<uint64_t>(sp.l_b()[ax] / sp.sigma());
//   // }
// }

// void write_last_cfg(void) {
//   std::string fn;
//   std::ostringstream oStrStream;
//   oStrStream << "read_init_cfg_pos_vel_Nm_" << sp.Nm() << ".xyz";
//   fn = oStrStream.str();
//   std::ofstream last_cfg_pos_vel(fn, std::ios::trunc);
//   last_cfg_pos_vel.precision(8);
//   for (int i = 0; i < sp.Nm(); ++i) {
//     last_cfg_pos_vel << r[i][0] << " " << r[i][1] << " " << r[i][2] << " "
//                      << v[i][0] << " " << v[i][1] << " " << v[i][2]
//                      << std::endl;
//   }
//   last_cfg_pos_vel.close();
// }

// void init_system(const int argc, const char **argv) {
//   // init_parameter();
//   // cell = new int64_t **[Cell_N[0] + 2];
//   // for (int64_t cx = 0; cx < Cell_N[0] + 2; cx++) {
//   //   cell[cx] = new int64_t *[Cell_N[1] + 2];
//   //   for (int64_t cy = 0; cy < Cell_N[1] + 2; cy++) {
//   //     cell[cx][cy] = new int64_t[Cell_N[2] + 2];
//   //   }
//   // }
//   // list = new int64_t[sp.Nm()];

//   // sp.read_arg(argc, argv);
//   std::cout << "Nm " << sp.Nm() << " Density " << sp.density() << " box "
//             << sp.l_b()[0] << std::endl;
//   // read_config();
//   std::cout << "Nm " << sp.Nm() << " Density " << sp.density() << " box "
//             << sp.l_b()[0] << std::endl;
//   init_position();
//   init_velocity();
//   init_force();
//   if (sp.gamma()) {
//     sp.calc_BD_factor();
//     init_gamma();
//   }
//   write_last_cfg();
// }

// void close_system(void) {}
