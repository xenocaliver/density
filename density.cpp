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
#include <cstdlib>
#include <cmath>
#include <thread>
#include <vector>
#include <boost/property_tree/json_parser.hpp>

extern "C" {
#include <fftw3.h>
}

#include "fourier_gauss.hpp"

namespace density {

extern fftw_complex* update_variable_node(std::vector<fftw_complex*>, fftw_complex*, uint64_t);
extern fftw_complex* update_check_node(std::vector<fftw_complex*>, uint64_t);

template<typename T> void fftw_memcpy(T* dest, T* src, uint64_t number_of_elements) {
    uint64_t uli;

    for(uli = 0; uli < number_of_elements; uli++) {
        dest[uli][0] = src[uli][0];
        dest[uli][1] = src[uli][1];
    }
}
}

int main(int argc, char* argv[]) {
    fftw_complex* ptr_fftw_complex;
    uint64_t uli;
    uint64_t vector_size = 1024;
    std::vector<fftw_complex*> variable_node_inputs(2, nullptr);
    std::vector<fftw_complex*> check_node_inputs(5, nullptr);
    fftw_complex* channel_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    double initial_dist = 1.0/(double)vector_size;

    for(uli = 0; uli < 2; uli++) {
        ptr_fftw_complex = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
        variable_node_inputs[uli] = ptr_fftw_complex;
    }
    for(uli = 0; uli < 5; uli++) {
        ptr_fftw_complex = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
        check_node_inputs[uli] = ptr_fftw_complex;
    }

    /* initial distribution */
    for(uli = 0; uli < vector_size; uli++) {
        variable_node_inputs[0][uli][0] = initial_dist; 
    }
    density::fftw_memcpy<fftw_complex>(variable_node_inputs[1], variable_node_inputs[0], vector_size);
    density::fftw_memcpy<fftw_complex>(variable_node_inputs[2], variable_node_inputs[0], vector_size);

    /* channel inputs */
    for(uli = 0; uli < vector_size; uli++) channel_input[uli][0] = density::fourier_gauss<double>((double)uli/vector_size, 0.5);

    /*variable node calculation */
    ptr_fftw_complex = density::update_variable_node(variable_node_inputs, channel_input, vector_size);
    
    /* copy probability distribution */
    for(uli = 0; uli < 5; uli++) density::fftw_memcpy<fftw_complex>(check_node_inputs[uli], ptr_fftw_complex, vector_size);

    /* free memrmoy */
    fftw_free(ptr_fftw_complex);

    /* check node calculation */
    ptr_fftw_complex = density::update_check_node(check_node_inputs, vector_size);

    fftw_free(ptr_fftw_complex);
    return(EXIT_SUCCESS);
}