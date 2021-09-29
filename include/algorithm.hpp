/**
 * @file algorithm.hpp
 * @author Haoran Chen (chen950302@live.com)
 * @brief
 * @version 0.1
 * @date 2021-09-18
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_
#include <cmath>  // floor
// #include <vector>
// #include <memory>
#include "init.hpp"

void cell_list(void);

double minium_image(const uint64_t &i, const uint64_t &j, const int &ax);

double LJ(uint64_t i, uint64_t j);
void MD_Step(void);
void calc_force(void);
void vel_correcter(void);


#endif  // ALGORITHM_HPP_
