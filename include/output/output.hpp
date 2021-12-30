/**
 * @file output.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef OUTPUT_OUTPUT_HPP_
#define OUTPUT_OUTPUT_HPP_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "initialization.hpp"

class output : virtual public initialization {
 private:
  std::string fn;
  std::ostringstream oss;
  int cfgs_pro_file = 1000;
  int cfgs_in_file = 0;
  std::ofstream last_cfg_pos_vel;
  std::ofstream cfg_data_file;
  std::ofstream cfg_time_file;
  std::ofstream energy_file;
  int cfg_file_array_Nm = 0;
  uint64_t cfg_freq_;
  uint64_t cfg_freq();

 public:
  output(/* args */);
  void print_energy();
  void write_cfg();
  void write_energy();
  void write_last_cfg();
  ~output();
};

// void print_Energy(void);

#endif  // OUTPUT_OUTPUT_HPP_
