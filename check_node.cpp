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
#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>

extern "C" {
#include <complex.h>
#include <fftw3.h>
}

namespace density {
extern double* get_two_way_tanh_map(double*, double*, uint64_t);
}

namespace density {

double* update_check_node(std::vector<double*> input, uint64_t vector_size) {
    uint64_t degree = input.size();
    double* ptr_output;
    double* ptr_tmp;
    uint64_t uli;

    ptr_output = get_two_way_tanh_map(input[0], input[1], vector_size);
    for(uli = 2; uli < degree; uli++) {
        ptr_tmp = get_two_way_tanh_map(input[uli], ptr_output, vector_size);
        fftw_free(ptr_output);
        ptr_output = ptr_tmp;
    }
    return(ptr_output);
}
}