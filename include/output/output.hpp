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

#include <string>
#include <iostream>

#include "initialization.hpp"

class output : virtual public initialization{
 private:
  std::string fn;
  std::ostringstream oStrStream;
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
