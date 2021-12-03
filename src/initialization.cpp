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

#include "initialization.hpp"

std::vector<std::array<double, 3>> r;
std::vector<std::array<double, 3>> dr;
std::vector<std::array<double, 3>> v;
std::vector<std::array<double, 3>> f0;
std::vector<std::array<double, 3>> f1;
std::vector<std::array<double, 3>> g0;
std::vector<std::array<double, 3>> g1;

void sys_param::read_arg(const int argc, const char **argv) {
  std::cout << "there" << std::endl;
  MD_time_ = std::stod(argv[1]);
  h_ = std::stod(argv[2]);
  density_ = std::stod(argv[3]);
  gamma_ = std::stod(argv[4]);

  MD_Steps_ = static_cast<uint64_t>(MD_time_ / h_);
  half_h_ = 0.5 * h_;
  half_h2_ = 0.5 * h_ * h_;
  time_100_ = static_cast<uint64_t>(100. / h_);
  time_10_ = static_cast<uint64_t>(10. / h_);
  time_1_ = static_cast<uint64_t>(1. / h_);
  time_01_ = static_cast<uint64_t>(0.1 / h_);
  time_001_ = static_cast<uint64_t>(0.01 / h_);
  time_0001_ = static_cast<uint64_t>(0.001 / h_);
}

uint64_t sys_param::MD_Steps() const { return MD_Steps_; }
void sys_param::Relax_Steps(const uint64_t input) { Relax_Steps_ = input; }
uint64_t sys_param::Relax_Steps() const { return Relax_Steps_; }

/**
 * @brief out put time step h
 *
 * @return double
 */
double sys_param::h() const { return h_; }
/**
 * @brief out put half time step
 *
 * @return double
 */
double sys_param::half_h() const { return half_h_; }
double sys_param::half_h2() const { return half_h2_; }
uint64_t sys_param::time_0001() const { return time_0001_; }
uint64_t sys_param::time_001() const { return time_001_; }
uint64_t sys_param::time_01() const { return time_01_; }
uint64_t sys_param::time_1() const { return time_1_; }
uint64_t sys_param::time_10() const { return time_10_; }
uint64_t sys_param::time_100() const { return time_100_; }

void sys_param::kT(const double input) { kT_ = input; }
double sys_param::kT() const { return kT_; }
void sys_param::m(const double input) { m_ = input; }
double sys_param::m() const { return m_; }

/**
 * @brief set the periodic boundary condition
 *
 * @param l_b_x
 * @param l_b_y
 * @param l_b_z
 */
void sys_param::l_b(const int ax, const double input) {
  l_b_[ax] = input;
  half_l_b_[ax] = 0.5 * input;
  inv_l_b_[ax] = 1. / input;
}
std::array<double, 3> sys_param::l_b() const { return l_b_; }
std::array<double, 3> sys_param::half_l_b() const { return half_l_b_; }
std::array<double, 3> sys_param::inv_l_b() const { return inv_l_b_; }

void sys_param::sigma(const double input) {
  sigma_ = input;
  sig2_ = input * input;
  r2_cut_ = input * input * pow(2., 1. / 3.);
}
double sys_param::sigma() const { return sigma_; }
double sys_param::r2_cut() const { return r2_cut_; }
double sys_param::sig2() const { return sig2_; }

void sys_param::epsilon(const double input) {
  epsilon_ = input;
}  // TODO e = kT
double sys_param::epsilon() const { return epsilon_; }

void sys_param::Nm(const uint64_t input) {
  Nm_ = input;
  density_ =
      static_cast<double>(input) * M_PI_4 * sig2_ / l_b_[0] / l_b_[1] / l_b_[2];
}

void sys_param::calc_Nm() {
  Nm_ = static_cast<uint64_t>(l_b_[0] * l_b_[1] * l_b_[2] * density_ / M_PI_4 /
                              sig2_);
}
uint64_t sys_param::Nm() const { return Nm_; }
double sys_param::density() const { return density_; }

double sys_param::C(const double x) {
  return 2. * x - 3. + 4. * exp(-x) - exp(-2. * x);
}
double sys_param::G(const double x) { return exp(x) - 2. * x - exp(-x); }

void sys_param::calc_BD_factor() {
  double gh = gamma_ * h_;

  double C_gh = C(gh);
  double G_gh = G(gh);
  double E_gh = -C(gh) * C(-gh) - pow(G(gh), 2);

  BD_g0_1_ = sqrt(kT_ / m_ / gamma_ / gamma_ * C_gh);
  BD_g1_1_ = sqrt(kT_ / m_ / gamma_ / gamma_ * E_gh / C_gh);
  std::cout << "E_gh " << E_gh << std::endl;
  BD_g1_2_ = G_gh / C_gh;

  BD_r_1_ = (1. - exp(-gh)) / gamma_;
  BD_r_2_ = (gh - 1. + exp(-gh)) / gamma_ / gamma_;

  BD_v_1_ = gamma_ / (exp(gh) - 1.);
  BD_v_2_ = (gh - 1. + exp(-gh)) / gamma_ / (exp(gh) - 1.);
  BD_v_3_ = (-gh - 1. + exp(gh)) / gamma_ / (exp(gh) - 1.);
}
double sys_param::gamma() const { return gamma_; }
double sys_param::BD_r_1() const { return BD_r_1_; }
double sys_param::BD_r_2() const { return BD_r_2_; }
double sys_param::BD_v_1() const { return BD_v_1_; }
double sys_param::BD_v_2() const { return BD_v_2_; }
double sys_param::BD_v_3() const { return BD_v_3_; }
double sys_param::BD_g0_1() const { return BD_g0_1_; }
double sys_param::BD_g1_1() const { return BD_g1_1_; }
double sys_param::BD_g1_2() const { return BD_g1_2_; }
sys_param sp;

int64_t step;

// TODO change the interpretation of Cell list
// uint64_t set_Cell_N[3];
// const uint64_t *const Cell_N = set_Cell_N;

// double set_Cell_l[3];
// const double *const Cell_l = set_Cell_l;

// int64_t ***cell;
// int64_t *list;

double E_kin, E_pot;
int print_E;
int FREQ_print_E = 10;

void read_config() {
  std::cout << "read initial file:" << std::endl;

  std::ifstream config_file;
  std::string line;
  // std::string head;
  int split_at;
  config_file.open("config.txt", std::ios::in);
  if (config_file.is_open()) {
    while (std::getline(config_file, line)) {
      if (line[0] != '#') {
        split_at = line.find('=');
        std::cout << line.substr(0, split_at - 1) << ":\t"
                  << line.substr(split_at + 2) << std::endl;

        const std::string &head = line.substr(0, split_at - 1);
        const std::string &value = line.substr(split_at + 2);

        if (head == "kT") {
          sp.kT(std::stod(value));
        } else if (head == "lx") {
          sp.l_b(0, std::stod(value));
        } else if (head == "ly") {
          sp.l_b(1, std::stod(value));
        } else if (head == "lz") {
          sp.l_b(2, std::stod(value));
        } else if (head == "sigma") {
          sp.sigma(std::stod(value));
        } else if (head == "epsilon") {
          sp.epsilon(std::stod(value));
        } else if (head == "m") {
          sp.m(std::stod(value));
        }
      }
    }
    sp.calc_Nm();
  }
}

void init_parameter(void) {
  sp.calc_BD_factor();
  // for (int ax = 0; ax < 3; ax++) {
  //   set_Cell_N[ax] = static_cast<uint64_t>(sp.l_b()[ax] / sp.sigma());
  // }
}

void init_position(void) {
  int row_x, row_y, row_z;
  row_x = static_cast<int>(sp.l_b()[0] / sp.sigma());
  row_y = static_cast<int>(sp.l_b()[1] / sp.sigma());
  row_z = static_cast<int>(sp.l_b()[2] / sp.sigma());
  int i = 0;
  std::array<double, 3> new_r;

  for (int r_z = 0; r_z < row_z; r_z++) {
    for (int r_y = 0; r_y < row_y; r_y++) {
      for (int r_x = 0; r_x < row_x; r_x++) {
        new_r[0] = r_x * sp.sigma();
        new_r[1] = r_y * sp.sigma();
        new_r[2] = r_z * sp.sigma();

        r.push_back(new_r);
        dr.push_back({0, 0, 0});
        if (++i >= sp.Nm()) {
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

  std::array<double, 3> new_g0, new_g1;

  for (uint64_t i = 0; i < sp.Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      new_g0[ax] = sp.BD_g0_1() * n_d(gen);
      new_g1[ax] = sp.BD_g1_1() * n_d(gen) + sp.BD_g1_2() * new_g0[ax];
    }
    g0.push_back(new_g0);
    g1.push_back(new_g1);
  }
}

void init_velocity(void) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, sp.kT());
  uint64_t i;
  int ax, th;
  double v_sum[3] = {0, 0, 0};

  std::array<double, 3> new_v;
  for (i = 0; i < sp.Nm(); i++) {
    new_v[0] = n_d(gen);
    new_v[1] = n_d(gen);
    new_v[2] = n_d(gen);
    v.push_back(new_v);
  }
#pragma omp parallel for
  for (i = 0; i < sp.Nm(); i++) {
    for (ax = 0; ax < 3; ax++) {
      v_sum[ax] += v[i][ax];
    }
  }
  std::cout << "initialization of velocity finished" << std::endl;
}

void init_force(void) {
  std::array<double, 3> new_f{0, 0, 0};
  for (uint64_t i = 0; i < sp.Nm(); i++) {
    f0.push_back(new_f);
    f1.push_back(new_f);
  }
}

void write_last_cfg(void) {
  std::string fn;
  std::ostringstream oStrStream;
  oStrStream << "read_init_cfg_pos_vel_Nm_" << sp.Nm() << ".xyz";
  fn = oStrStream.str();
  std::ofstream last_cfg_pos_vel(fn, std::ios::trunc);
  last_cfg_pos_vel.precision(8);
  for (int i = 0; i < sp.Nm(); ++i) {
    last_cfg_pos_vel << r[i][0] << " " << r[i][1] << " " << r[i][2] << " "
                     << v[i][0] << " " << v[i][1] << " " << v[i][2]
                     << std::endl;
  }
  last_cfg_pos_vel.close();
}

void init_system(const int argc, const char **argv) {
  // init_parameter();
  // cell = new int64_t **[Cell_N[0] + 2];
  // for (int64_t cx = 0; cx < Cell_N[0] + 2; cx++) {
  //   cell[cx] = new int64_t *[Cell_N[1] + 2];
  //   for (int64_t cy = 0; cy < Cell_N[1] + 2; cy++) {
  //     cell[cx][cy] = new int64_t[Cell_N[2] + 2];
  //   }
  // }
  // list = new int64_t[sp.Nm()];

  sp.read_arg(argc, argv);
  read_config();

  init_position();
  init_velocity();
  init_force();
  if (sp.gamma()) {
    sp.calc_BD_factor();
    init_gamma();
  }
  write_last_cfg();
}

void close_system(void) {}
