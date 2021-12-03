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

extern std::vector<std::array<double, 3>> r;
extern std::vector<std::array<double, 3>> dr;
extern std::vector<std::array<double, 3>> v;
extern std::vector<std::array<double, 3>> f0;
extern std::vector<std::array<double, 3>> f1;
extern std::vector<std::array<double, 3>> g0;
extern std::vector<std::array<double, 3>> g1;

class sys_param {
 private:
  double MD_time_;
  double Relax_time_;
  double h_;
  double half_h_;
  double half_h2_;
  uint64_t MD_Steps_;
  uint64_t Relax_Steps_;
  uint64_t step_;
  uint64_t time_0001_;
  uint64_t time_001_;
  uint64_t time_01_;
  uint64_t time_1_;
  uint64_t time_10_;
  uint64_t time_100_;

  double kT_;
  double m_;
  double density_;
  std::array<double, 3> l_b_;
  std::array<double, 3> half_l_b_;
  std::array<double, 3> inv_l_b_;
  double Nm_;
  double gamma_;
  double sigma_;
  double epsilon_;
  double r2_cut_;
  double sig2_;

  double BD_r_1_;
  double BD_r_2_;
  double BD_v_1_;
  double BD_v_2_;
  double BD_v_3_;
  double BD_g0_1_;
  double BD_g1_1_;
  double BD_g1_2_;
  double C(const double x);
  double G(const double x);

 public:
  void read_arg(const int argc, const char **argv);

  void MD_Steps(const double input);
  uint64_t MD_Steps() const;
  void Relax_Steps(const uint64_t input);
  uint64_t Relax_Steps() const;
  void step(const double input);
  uint64_t step() const;
  void MD_time(const double input);
  double MD_time() const;

  void kT(const double input);
  double kT() const;
  void m(const double input);
  double m() const;

  void h(const double input);
  double h() const;
  double half_h() const;
  double half_h2() const;
  uint64_t time_0001() const;
  uint64_t time_001() const;
  uint64_t time_01() const;
  uint64_t time_1() const;
  uint64_t time_10() const;
  uint64_t time_100() const;

  void l_b(const int ax, const double input);
  std::array<double, 3> l_b() const;
  std::array<double, 3> half_l_b() const;
  std::array<double, 3> inv_l_b() const;

  void Nm(const uint64_t input);
  uint64_t Nm() const;
  void calc_Nm();
  double density() const;

  void gamma(const double input);
  double gamma() const;
  void sigma(const double input);
  double sigma() const;
  void epsilon(const double input);
  double epsilon() const;

  double r2_cut() const;
  double sig2() const;

  void calc_BD_factor();
  double BD_r_1() const;
  double BD_r_2() const;
  double BD_v_1() const;
  double BD_v_2() const;
  double BD_v_3() const;
  double BD_g0_1() const;
  double BD_g1_1() const;
  double BD_g1_2() const;

  double random() const;
};

extern sys_param sp;

extern int64_t step;

// extern uint64_t set_Cell_N[3];
// extern const uint64_t *const Cell_N;

// extern double set_Cell_l[3];
// extern const double *const Cell_l;

// extern int64_t ***cell;
// extern int64_t *list;

extern double E_kin, E_pot;
extern int print_E;
extern int FREQ_print_E;

void read_config(void);
void init_parameter(void);
void init_gamma(void);
void init_position(void);
void init_velocity(void);
void init_system(const int argc, const char **argv);
void close_system(void);
void write_last_cfg(void);

#endif  // INITIALIZATION_HPP_
