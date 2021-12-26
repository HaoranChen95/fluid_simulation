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
