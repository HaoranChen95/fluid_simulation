/**
 * @file init.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-16
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef INITIALIZATION_INITIALIZATION_HPP_
#define INITIALIZATION_INITIALIZATION_HPP_

#ifdef _OPENMP
#define N_THREADS 4

#else
#define N_THREADS 1
#endif  // _OPENMP

#include <omp.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "box.hpp"
#include "fluctuation.hpp"
#include "force.hpp"
#include "velocity.hpp"

class initialization : protected fluctuation,
                       protected velocity,
                       protected force {
 private:
 public:
  // initialization();
  void init(const int argc, const char **argv);
  void read_arg(const int argc, const char **argv);
  void read_config();
};

// void init_system(const int argc, const char **argv);
// void close_system(void);
// void write_last_cfg(void);

#endif  // INITIALIZATION_INITIALIZATION_HPP_
