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

#ifndef INIT_HPP_
#define INIT_HPP_

#ifdef _OPENMP
#define N_THREADS 4

#else
#define N_THREADS 1
#endif  // _OPENMP

#include <omp.h>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

extern double **r0;
extern double **r1;
extern double **v;
extern double **f0;
extern double **f1;
extern double ***MP_f;

extern double set_kT;
extern const double &kT;

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
extern double set_dt;
extern const double &dt;
extern double set_half_dt;
extern const double &half_dt;
extern double set_half_dt2;
extern const double &half_dt2;

extern uint64_t set_Cell_N[3];
extern const uint64_t *const Cell_N;

extern double set_Cell_l[3];
extern const double *const Cell_l;

extern int64_t ***cell;
extern int64_t *list;

extern double E_kin, E_pot;

extern double set_r2_cut;
extern const double &r2_cut;
extern double set_sig2;
extern const double &sig2;

void read_config(void);
void init_parameter(void);
void init_position(void);
void init_velocity(void);
void init_system(void);
void close_system(void);
void write_last_cfg(void);

#endif  // INIT_HPP_
