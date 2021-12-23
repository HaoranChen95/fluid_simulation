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

initialization::initialization(const int argc, const char **argv) {
  read_arg(argc, argv);
  read_config();
}

void initialization::read_arg(const int argc, const char **argv) {
  MD_time(std::stod(argv[1]));
  h(std::stod(argv[2]));
  density(std::stod(argv[3]));
  gamma(std::stod(argv[4]));
}

void initialization::read_config() {
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
          kT(std::stod(value));
        } else if (head == "lx") {
          l_b(0, std::stod(value));
        } else if (head == "ly") {
          l_b(1, std::stod(value));
        } else if (head == "lz") {
          l_b(2, std::stod(value));
        } else if (head == "sigma") {
          sigma(std::stod(value));
        } else if (head == "epsilon") {
          epsilon(std::stod(value));
        } else if (head == "m") {
          m(std::stod(value));
        }
      }
    }
    calc_Nm();
    calc_density();
  }
}

// int argc_;
// const char** argv_;
// initialization sp(argc_, argv_);

// int64_t step;

// uint64_t set_Cell_N[3];
// const uint64_t *const Cell_N = set_Cell_N;

// double set_Cell_l[3];
// const double *const Cell_l = set_Cell_l;

// int64_t ***cell;
// int64_t *list;

// double E_kin, E_pot;
// int print_E;
// int FREQ_print_E = 10;

// void init_parameter(void) {
//   sp.calc_BD_factor();
//   // for (int ax = 0; ax < 3; ax++) {
//   //   set_Cell_N[ax] = static_cast<uint64_t>(sp.l_b()[ax] / sp.sigma());
//   // }
// }

// void init_position(void) {
//   int row_x, row_y, row_z;
//   row_x = static_cast<int>(sp.l_b()[0] / sp.sigma());
//   row_y = static_cast<int>(sp.l_b()[1] / sp.sigma());
//   row_z = static_cast<int>(sp.l_b()[2] / sp.sigma());
//   int i = 0;
//   std::array<double, 3> new_r;

//   for (int r_z = 0; r_z < row_z; r_z++) {
//     for (int r_y = 0; r_y < row_y; r_y++) {
//       for (int r_x = 0; r_x < row_x; r_x++) {
//         new_r[0] = r_x * sp.sigma();
//         new_r[1] = r_y * sp.sigma();
//         new_r[2] = r_z * sp.sigma();

//         r.push_back(new_r);
//         dr.push_back({0, 0, 0});
//         if (++i >= sp.Nm()) {
//           goto finish;
//         }
//       }
//     }
//   }
// finish:
//   std::cout << "initialization of position finished" << std::endl;
// }

// void init_gamma(void) {
//   std::random_device rd{};
//   std::mt19937 gen{rd()};
//   std::normal_distribution<double> n_d(0.0, 1.);

//   std::array<double, 3> new_g0, new_g1;

//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     for (int ax = 0; ax < 3; ax++) {
//       new_g0[ax] = sp.BD_g0_1() * n_d(gen);
//       new_g1[ax] = sp.BD_g1_1() * n_d(gen) + sp.BD_g1_2() * new_g0[ax];
//     }
//     g0.push_back(new_g0);
//     g1.push_back(new_g1);
//   }
// }

// void init_velocity(void) {
//   std::random_device rd{};
//   std::mt19937 gen{rd()};
//   std::normal_distribution<double> n_d(0.0, sp.kT());
//   uint64_t i;
//   int ax, th;
//   double v_sum[3] = {0, 0, 0};

//   std::array<double, 3> new_v;
//   for (i = 0; i < sp.Nm(); i++) {
//     new_v[0] = n_d(gen);
//     new_v[1] = n_d(gen);
//     new_v[2] = n_d(gen);
//     v.push_back(new_v);
//   }
// #pragma omp parallel for
//   for (i = 0; i < sp.Nm(); i++) {
//     for (ax = 0; ax < 3; ax++) {
//       v_sum[ax] += v[i][ax];
//     }
//   }
//   std::cout << "initialization of velocity finished" << std::endl;
// }

// void init_force(void) {
//   std::array<double, 3> new_f{0, 0, 0};
//   for (uint64_t i = 0; i < sp.Nm(); i++) {
//     f0.push_back(new_f);
//     f1.push_back(new_f);
//   }
// }

// void write_last_cfg(void) {
//   std::string fn;
//   std::ostringstream oStrStream;
//   oStrStream << "read_init_cfg_pos_vel_Nm_" << sp.Nm() << ".xyz";
//   fn = oStrStream.str();
//   std::ofstream last_cfg_pos_vel(fn, std::ios::trunc);
//   last_cfg_pos_vel.precision(8);
//   for (int i = 0; i < sp.Nm(); ++i) {
//     last_cfg_pos_vel << r[i][0] << " " << r[i][1] << " " << r[i][2] << " "
//                      << v[i][0] << " " << v[i][1] << " " << v[i][2]
//                      << std::endl;
//   }
//   last_cfg_pos_vel.close();
// }

// void init_system(const int argc, const char **argv) {
//   // init_parameter();
//   // cell = new int64_t **[Cell_N[0] + 2];
//   // for (int64_t cx = 0; cx < Cell_N[0] + 2; cx++) {
//   //   cell[cx] = new int64_t *[Cell_N[1] + 2];
//   //   for (int64_t cy = 0; cy < Cell_N[1] + 2; cy++) {
//   //     cell[cx][cy] = new int64_t[Cell_N[2] + 2];
//   //   }
//   // }
//   // list = new int64_t[sp.Nm()];

//   // sp.read_arg(argc, argv);
//   std::cout << "Nm " << sp.Nm() << " Density " << sp.density() << " box "
//             << sp.l_b()[0] << std::endl;
//   // read_config();
//   std::cout << "Nm " << sp.Nm() << " Density " << sp.density() << " box "
//             << sp.l_b()[0] << std::endl;
//   init_position();
//   init_velocity();
//   init_force();
//   if (sp.gamma()) {
//     sp.calc_BD_factor();
//     init_gamma();
//   }
//   write_last_cfg();
// }

// void close_system(void) {}
