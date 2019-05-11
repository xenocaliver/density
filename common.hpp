/* 
 * This file is part of the density distribution (https://github.com/xenocaliver/density).
 * Copyright (c) 2019 Akiyoshi Hashimoto.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>
 */

#include <cstdlib>
#include <cstdint>

constexpr uint64_t vector_size = 16384;
constexpr int64_t upper_bound = (vector_size - 1)/2;
constexpr int64_t lower_bound = -(int64_t)(vector_size - 1)/2;
constexpr int64_t half_upper_bound = (vector_size - 1)/4;
constexpr int64_t half_lower_bound = -(int64_t)(vector_size - 1)/4;
constexpr double delta = 0.05;
constexpr double error_probability_threshold = 1.0e-6;