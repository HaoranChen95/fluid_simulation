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

#ifndef INITIALIZATION_HPP_
#define INITIALIZATION_HPP_

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
#include "brown_factor.hpp"
#include "particle_parameter.hpp"
#include "time_step.hpp"
#include "velocity.hpp"
#include "force.hpp"

extern std::vector<std::array<double, 3>> g0;
extern std::vector<std::array<double, 3>> g1;

class initialization : public brown_factor, protected velocity, protected force {
 private:
 public:
  initialization(const int argc, const char **argv);
  void read_arg(const int argc, const char **argv);
  void read_config();
};

// extern int argc_;
// extern const char** argv_;

// extern initialization sp;

// extern int64_t step;

// extern uint64_t set_Cell_N[3];
// extern const uint64_t *const Cell_N;

// extern double set_Cell_l[3];
// extern const double *const Cell_l;

// extern int64_t ***cell;
// extern int64_t *list;

extern double E_kin, E_pot;
extern int print_E;
extern int FREQ_print_E;

// void read_config(void);
// void init_parameter(void);
// void init_gamma(void);
// void init_position(void);
// void init_velocity(void);
// void init_system(const int argc, const char **argv);
// void close_system(void);
// void write_last_cfg(void);

#endif  // INITIALIZATION_HPP_
