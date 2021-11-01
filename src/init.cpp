/**
 * @file init.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-16
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "init.hpp"

double **dr; /** @brief difference of position */
double **r1; /** @brief difference of position */
double **v;  /** @brief velocity */
double **f0; /** @brief previos force */
double **f1; /** @brief new force */
double **g0; /** @brief previos random number */
double **g1; /** @brief new random number */

double set_kT;
const double &kT = set_kT;

double set_kT_2gamma_over_m;
double &kT_2gamma_over_m = set_kT_2gamma_over_m;

double set_m;
double &m = set_m;

bool open_fluid;

double set_l_b[3];
const double *const l_b = set_l_b;

double set_half_l_b[3];
const double *const half_l_b = set_half_l_b;

double set_inv_l_b[3];
const double *const inv_l_b = set_inv_l_b;

uint64_t set_Nm;
const uint64_t &Nm = set_Nm;
double set_gamma;
double set_sigma;
double set_epsilon;
double set_density;
const double &gam = set_gamma;
const double &sig = set_sigma;
const double &eps = set_epsilon;
const double &density = set_density;

double set_r2_cut;
const double &r2_cut = set_r2_cut;
double set_sig2;
const double &sig2 = set_sig2;

int64_t MD_Steps;
int64_t Relax_Steps;
int64_t step;
double MD_time;

uint64_t set_time_1;
const uint64_t &time_1 = set_time_1;
uint64_t set_time_01;
const uint64_t &time_01 = set_time_01;
uint64_t set_time_001;
const uint64_t &time_001 = set_time_001;
uint64_t set_time_0001;
const uint64_t &time_0001 = set_time_0001;
uint64_t set_time_10;
const uint64_t &time_10 = set_time_10;
uint64_t set_time_100;
const uint64_t &time_100 = set_time_100;
uint64_t set_time_1000;
const uint64_t &time_1000 = set_time_1000;

double set_dt;
const double &dt = set_dt;

double set_half_dt;
const double &half_dt = set_half_dt;

double set_half_dt2;
const double &half_dt2 = set_half_dt2;

uint64_t set_Cell_N[3];
const uint64_t *const Cell_N = set_Cell_N;

double set_Cell_l[3];
const double *const Cell_l = set_Cell_l;

int64_t ***cell;
int64_t *list;

double E_kin, E_pot;
int print_E;
int FREQ_print_E = 10;

// randem constant part
double C_gamh;
double G_gamh;
double E_gamh;
double set_const_g0_1;
double set_const_g1_1;
double set_const_g1_2;
const double &const_g0_1 = set_const_g0_1;
const double &const_g1_1 = set_const_g1_1;
const double &const_g1_2 = set_const_g1_2;
double set_const_r_1;
double set_const_r_2;
double set_const_v_1;
double set_const_v_2;
double set_const_v_3;
const double &const_r_1 = set_const_r_1;
const double &const_r_2 = set_const_r_2;
const double &const_v_1 = set_const_v_1;
const double &const_v_2 = set_const_v_2;
const double &const_v_3 = set_const_v_3;

void read_config() {
  std::cout << "read initial file:" << std::endl;
  std::cout << MD_Steps << std::endl;

  std::ifstream config_file;
  std::string line;
  // std::string head;
  int split_at;
  config_file.open("config.txt", std::ios::in);
  if (config_file.is_open()) {
    while (std::getline(config_file, line)) {
      if (line[0] != '#') {
        split_at = line.find('=');
        std::cout << line.substr(0, split_at - 1) << "-"
                  << line.substr(split_at + 2) << std::endl;

        const std::string &head = line.substr(0, split_at - 1);
        const std::string &value = line.substr(split_at + 2);

        if (head == "kT") {
          set_kT = std::stod(value);
        } else if (head == "gamma") {
          set_gamma = std::stod(value);
        } else if (head == "lx") {
          set_l_b[0] = std::stod(value);
        } else if (head == "ly") {
          set_l_b[1] = std::stod(value);
        } else if (head == "lz") {
          set_l_b[2] = std::stod(value);
        } else if (head == "density") {
          set_density = std::stod(value);
        } else if (head == "sigma") {
          set_sigma = std::stod(value);
        } else if (head == "epsilon") {
          set_epsilon = std::stod(value);
        } else if (head == "m") {
          set_m = std::stod(value);
        } else if (head == "open_fluid") {
          open_fluid = std::stoi(value);
        }
        // else if (head == "step length") {
        //   set_dt = std::stod(value);
        // }
      }
    }
  }
}

void init_parameter(void) {
  set_Nm = static_cast<uint64_t>(density * l_b[0] * l_b[1] * l_b[2]);

  for (int ax = 0; ax < 3; ax++) {
    set_half_l_b[ax] = 0.5 * l_b[ax];
    set_inv_l_b[ax] = 1. / l_b[ax];
    set_Cell_N[ax] = static_cast<uint64_t>(l_b[ax] / sig);
  }

  set_r2_cut = sig * sig * pow(2., 1. / 3.);
  set_sig2 = sig * sig;
  set_half_dt = 0.5 * dt;
  set_half_dt2 = 0.5 * dt * dt;

  set_kT_2gamma_over_m = 2. * gam * kT / m;
  double gamh = gam * dt;
  C_gamh = 2. * gamh - 3. + 4. * exp(-gamh) - exp(-2. * gamh);
  G_gamh = exp(gamh) - 2. * gamh - exp(-gamh);
  E_gamh = 16. * (exp(gamh) + exp(-gamh)) -
           4. * (exp(2. * gamh) + exp(-2. * gamh)) -
           4. * gamh * (exp(gamh) - exp(-gamh)) +
           2. * gamh * (exp(2. * gamh) - exp(-2. * gamh)) - 24.;
  set_const_g0_1 = sqrt(kT / m / gam / gam * C_gamh);
  set_const_g1_1 = sqrt(kT / m / gam / gam * E_gamh / C_gamh);
  set_const_g1_2 = G_gamh / C_gamh;
  set_const_r_1 = (1. - exp(-gamh)) / gam;
  set_const_r_2 = (gamh - 1. + exp(-gamh)) / gam / gam;
  set_const_v_1 = gam / (exp(gamh) - 1.);
  set_const_v_2 = (gamh - 1. + exp(-gamh)) / gam / (exp(gamh) - 1.);
  set_const_v_3 = (-gamh - 1. + exp(gamh)) / gam / (exp(gamh) - 1.);

  set_time_1000 = static_cast<uint64_t>(1000. / dt);
  set_time_100 = static_cast<uint64_t>(100. / dt);
  set_time_10 = static_cast<uint64_t>(10. / dt);
  set_time_1 = static_cast<uint64_t>(1. / dt);
  set_time_01 = static_cast<uint64_t>(0.1 / dt);
  set_time_001 = static_cast<uint64_t>(0.01 / dt);
  set_time_0001 = static_cast<uint64_t>(0.001 / dt);
}

void init_position(void) {
  int row_x, row_y, row_z;
  row_x = static_cast<int>(l_b[0] / sig);
  row_y = static_cast<int>(l_b[1] / sig);
  row_z = static_cast<int>(l_b[2] / sig);
  int i = 0;

  for (int r_z = 0; r_z < row_z; r_z++) {
    for (int r_y = 0; r_y < row_y; r_y++) {
      for (int r_x = 0; r_x < row_x; r_x++) {
        r1[0][i] = r_x * sig;
        r1[1][i] = r_y * sig;
        r1[2][i] = r_z * sig;
        if (++i >= Nm) {
          goto finish;
        }
      }
    }
  }
finish:
  std::cout << "initialization of position finished" << std::endl;
}

void init_gamma(void) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, 1.);

  for (uint64_t i = 0; i < Nm; i++) {
    for (int ax = 0; ax < 3; ax++) {
      g0[ax][i] = const_g0_1 * n_d(gen);
      g1[ax][i] = const_g1_1 * n_d(gen) + const_g1_2 * g0[ax][i];
    }
  }
}

void init_velocity(void) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, kT);
  uint64_t i;
  int ax, th;
  double v_sum[3] = {0, 0, 0};

  for (i = 0; i < Nm; i++) {
    v[0][i] = n_d(gen);
    v[1][i] = n_d(gen);
    v[2][i] = n_d(gen);
  }
#pragma omp parallel for
  for (i = 0; i < Nm; i++) {
    for (ax = 0; ax < 3; ax++) {
      v_sum[ax] += v[ax][i];
    }
  }

  // for (ax = 0; ax < 3; ax++) {
  //   v_sum[ax] /= static_cast<double>(Nm);
  // }

  // for (i = 0; i < Nm; i++) {
  //   for (ax = 0; ax < 3; ax++) {
  //     v[ax][i] -= v_sum[ax];
  //   }
  // }

  std::cout << "initialization of velocity finished" << std::endl;
}

void write_last_cfg(void) {
  std::string fn;
  std::ostringstream oStrStream;
  oStrStream << "read_init_cfg_pos_vel_Nm_" << Nm << ".xyz";
  fn = oStrStream.str();
  std::ofstream last_cfg_pos_vel(fn, std::ios::trunc);
  last_cfg_pos_vel.precision(8);
  for (int i = 0; i < Nm; ++i) {
    last_cfg_pos_vel << r1[0][i] << " " << r1[1][i] << " " << r1[2][i] << " "
                     << v[0][i] << " " << v[1][i] << " " << v[2][i]
                     << std::endl;
    //  << " "
    //  << ex[i] << " " << ey[i] << " " << ez[i] << endl;
  }
  last_cfg_pos_vel.close();
}

void init_system(void) {
  init_parameter();
  // set_Nm = 2;
  dr = new double *[3];
  r1 = new double *[3];
  v = new double *[3];
  f0 = new double *[3];
  f1 = new double *[3];
  g0 = new double *[3];
  g1 = new double *[3];

  for (int ax = 0; ax < 3; ax++) {
    dr[ax] = new double[Nm];
    r1[ax] = new double[Nm];
    v[ax] = new double[Nm];
    f0[ax] = new double[Nm];
    f1[ax] = new double[Nm];
    g0[ax] = new double[Nm];
    g1[ax] = new double[Nm];
  }

  cell = new int64_t **[Cell_N[0] + 2];
  for (int64_t cx = 0; cx < Cell_N[0] + 2; cx++) {
    cell[cx] = new int64_t *[Cell_N[1] + 2];
    for (int64_t cy = 0; cy < Cell_N[1] + 2; cy++) {
      cell[cx][cy] = new int64_t[Cell_N[2] + 2];
    }
  }

  list = new int64_t[Nm];
  init_position();
  init_velocity();
  write_last_cfg();
}

void close_system(void) {
  for (int ax = 0; ax < 3; ax++) {
    delete[] dr[ax];
    delete[] r1[ax];
    delete[] v[ax];
    delete[] f0[ax];
    delete[] f1[ax];
    delete[] g0[ax];
    delete[] g1[ax];
  }
  delete[] dr;
  delete[] r1;
  delete[] v;
  delete[] f0;
  delete[] f1;
  delete[] g0;
  delete[] g1;
}
