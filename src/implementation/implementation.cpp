/**
 * @file implementation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "implementation.hpp"

implementation::implementation() {}

void implementation::relaxation() {
  std::cout << "finished init" << std::endl;
  if (gamma()) {
    std::cout << "BD" << std::endl;
    run_BD_step();
  } else {
    std::cout << "MD" << std::endl;
    run_MD_step();
  }
}

implementation::~implementation() {}
