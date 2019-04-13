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
#include <cstdint>
#include <vector>
#include <cmath>

extern "C" {
#include <fftw3.h>
}

namespace density {

uint64_t quasi_gallager(uint64_t x, uint64_t y, uint64_t vector_size) {
    double a, b;
    double delta = 1.0/(double)vector_size;
    double z = 0.0;
    uint64_t rtnv;

    a = 0.5*((double)x)*delta;
    b = 0.5*((double)y)*delta;

    z = 2.0*atanh(tanh(a)*tanh(b));
    rtnv = (uint64_t)std::round(z/delta);
    return(rtnv);
}

fftw_complex* get_two_way_tanh_map(fftw_complex* a, fftw_complex* b, uint64_t vector_size) {
    uint64_t uli, ulj, ulk;
    fftw_complex* ptr_output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);

    for(uli = 0; uli < vector_size; uli++) {
        ptr_output[uli][0] = 0.0;
        ptr_output[uli][1] = 0.0;
    }

    for(uli = 0; uli < vector_size; uli++) {
        for(ulj = 0; ulj < vector_size; ulj++) {
            ulk = quasi_gallager(uli, ulj, vector_size);
            ptr_output[ulk][0] += a[uli][0]*b[uli][0] - a[uli][1]*b[uli][1];
            ptr_output[ulk][1] += a[uli][1]*b[uli][0] + a[uli][0]*b[uli][1];
        }
    }
    return(ptr_output);
}

fftw_complex* update_check_node(std::vector<fftw_complex*> input, uint64_t vector_size) {
    uint64_t degree = input.size();
    fftw_complex* ptr_output;
    fftw_complex* ptr_tmp;
    uint64_t uli;

    ptr_output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    for(uli = 0; uli < vector_size; uli++) {
        ptr_output[uli][0] = (double)0.0;
        ptr_output[uli][1] = (double)0.0;
    }

    ptr_output = get_two_way_tanh_map(input[0], input[1], vector_size);
    for(uli = 2; uli < degree; uli++) {
        ptr_tmp = get_two_way_tanh_map(input[uli], ptr_output, vector_size);
        fftw_free(ptr_output);
        ptr_output = ptr_tmp;
    }
    return(ptr_output);
}
}