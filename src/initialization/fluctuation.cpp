/**
 * @file fluctuation.cpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-12-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "fluctuation.hpp"

void fluctuation::init_fluctuation() {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, 1.);

  std::array<double, 3> g{0, 0, 0};

  for (uint64_t i = 0; i < Nm(); i++) {
    g0.push_back(g);
    g1.push_back(g);
  }
  generate_Gamma();
  std::cout << "initialization of fluctuation finished: size " << g0.size()
            << ", " << g1.size() << std::endl;
}

void fluctuation::generate_Gamma(void) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> n_d(0.0, 1.);

  for (uint64_t i = 0; i < Nm(); i++) {
    for (int ax = 0; ax < 3; ax++) {
      g0[i][ax] = BD_g0_1() * n_d(gen);
      g1[i][ax] = BD_g1_1() * n_d(gen) + BD_g1_2() * g0[i][ax];
    }
  }
}
