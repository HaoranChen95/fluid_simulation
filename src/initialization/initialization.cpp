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

void initialization::print_init() {
  print_time();
  print_particle();
  print_box();
  if (gamma()) {
    print_brown_factor();
  }
}

void initialization::initialize(const int argc, const char **argv) {
  read_arg(argc, argv);
  read_config();
  init_position();
  init_cell_list(sqrt(r2_cut()));
  init_velocity();
  init_force();
  if (gamma()) {
    calc_BD_factor();
    init_fluctuation();
  }
  print_init();
}

void initialization::read_arg(const int argc, const char **argv) {
  MD_time(std::stod(argv[1]));
  h(std::stod(argv[2]));
  density(std::stod(argv[3]));
  kT(std::stod(argv[4]));
  gamma(std::stod(argv[5]));
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

        if (head == "lx") {
          l_b(0, std::stod(value));
        } else if (head == "ly") {
          l_b(1, std::stod(value));
        } else if (head == "lz") {
          l_b(2, std::stod(value));
        } else if (head == "sigma") {
          sigma(std::stod(value));
        } else if (head == "m") {
          m(std::stod(value));
        } else if (head == "relax_time") {
          Relax_time(std::stod(value));
        }
      }
    }
    calc_Nm();
    calc_density();
  }
}
