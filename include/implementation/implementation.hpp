/**
 * @file implementation.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-26
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef IMPLEMENTATION_IMPLEMENTATION_HPP_
#define IMPLEMENTATION_IMPLEMENTATION_HPP_
#include <cmath>  // floor
// #include <vector>
// #include <memory>
#include "MD_step.hpp"

// void cell_list(void);
class implementation : protected MD_step {
 private:
  /* data */
 public:
  implementation(const int argc, const char **argv);
  // void init(const int argc, const char **argv);
  void relaxation();
  ~implementation();
};

// void MD_Step(void);
// void vel_correcter(void);

#endif  // IMPLEMENTATION_IMPLEMENTATION_HPP_
