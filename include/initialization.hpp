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
  double kT_;
  double kT_2gamma_over_m_;
  double m_;
  std::array<double, 3> l_b_;
  std::array<double, 3> half_l_b_;
  std::array<double, 3> inv_l_b_;
  double Nm_;
  double gamma_;
  double sigma_;
  double epsilon_;
  double density_;
  double r2_cut_;
  double sig2_;

  uint64_t MD_Step_;
  uint64_t Relax_Steps_;
  uint64_t step_;
  double MD_time_;
  double h_;
  double half_h_;
  double half_h2_;
  uint64_t time_0001_;
  uint64_t time_001_;
  uint64_t time_01_;
  uint64_t time_1_;
  uint64_t time_10_;
  uint64_t time_100_;

  double C_gamh_;
  double G_gamh_;
  double E_gamh_;
  double const_g0_1_;
  double const_g1_1_;
  double const_g1_2_;
  double const_r_1_;
  double const_r_2_;
  double const_v_1_;
  double const_v_2_;
  double const_v_3_;

 public:

  /* finished */
  sys_param(/* args */);
  void kT(const double input);
  double kT() const;
  void kT_2gamma_over_m(const double input);
  double kT_2gamma_over_m() const;
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

  /*haven't finished */

  void Nm(const double input);
  double Nm() const;
  void gamma(const double input);
  double gamma() const;
  void sigma(const double input);
  double sigma() const;
  void epsilon(const double input);
  double epsilon() const;
  void density(const double input);
  double density() const;

  std::array<double, 3> l_b() const;
  std::array<double, 3> half_l_b() const;
  std::array<double, 3> inv_l_b() const;
  void r2_cut(const double input);
  double r2_cut() const;
  void sig2(const double input);
  double sig2() const;
  void MD_Step(const double input);
  uint64_t MD_Step() const;
  void Relax_Steps(const double input);
  uint64_t Relax_Steps() const;
  void step(const double input);
  uint64_t step() const;
  void MD_time(const double input);
  double MD_time() const;
  void C_gamh(const double input);
  double C_gamh() const;
  void G_gamh(const double input);
  double G_gamh() const;
  void E_gamh(const double input);
  double E_gamh() const;
  void const_g0_1(const double input);
  double const_g0_1() const;
  void const_g1_1(const double input);
  double const_g1_1() const;
  void const_g1_2(const double input);
  double const_g1_2() const;
  void const_r_1(const double input);
  double const_r_1() const;
  void const_r_2(const double input);
  double const_r_2() const;
  void const_v_1(const double input);
  double const_v_1() const;
  void const_v_2(const double input);
  double const_v_2() const;
  void const_v_3(const double input);
  double const_v_3() const;
  ~sys_param();
};

extern sys_param sp;

extern double set_kT_2gamma_over_m;
extern double &kT_2gamma_over_m;

extern bool open_fluid;

extern double set_l_b[3];
extern const double *const l_b;

extern double set_half_l_b[3];
extern const double *const half_l_b;

extern double set_inv_l_b[3];
extern const double *const inv_l_b;

extern uint64_t set_Nm;
extern const uint64_t &Nm;
extern double set_gamma;
extern double set_sigma;
extern double set_epsilon;
extern double set_density;
extern const double &gam;
extern const double &sig;
extern const double &eps;
extern const double &density;

extern int64_t MD_Steps;
extern int64_t Relax_Steps;
extern int64_t step;
extern double MD_time;

extern uint64_t set_Cell_N[3];
extern const uint64_t *const Cell_N;

extern double set_Cell_l[3];
extern const double *const Cell_l;

extern int64_t ***cell;
extern int64_t *list;

extern double E_kin, E_pot;
extern int print_E;
extern int FREQ_print_E;

extern double set_r2_cut;
extern const double &r2_cut;
extern double set_sig2;
extern const double &sig2;

extern double set_const_r_1;
extern double set_const_r_2;
extern double set_const_v_1;
extern double set_const_v_2;
extern double set_const_v_3;
extern const double &const_r_1;
extern const double &const_r_2;
extern const double &const_v_1;
extern const double &const_v_2;
extern const double &const_v_3;

extern double set_const_g0_1;
extern double set_const_g1_1;
extern double set_const_g1_2;
extern const double &const_g0_1;
extern const double &const_g1_1;
extern const double &const_g1_2;

void read_arg(const int argc, const char **argv);
void read_config(void);
void init_parameter(void);
void init_gamma(void);
void init_position(void);
void init_velocity(void);
void init_system(void);
void close_system(void);
void write_last_cfg(void);

#endif  // INITIALIZATION_HPP_
