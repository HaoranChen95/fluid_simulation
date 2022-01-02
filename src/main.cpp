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
  fluid_simulation fluid(argc, argv);
  
  fluid.relax();
  fluid.implement();
  return 0;
}
