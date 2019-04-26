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
#include <cstring>
#include <complex>

extern "C" {
#include <fftw3.h>
}

#include "transform_complex.hpp"

namespace density {

double* update_variable_node(std::vector<fftw_complex*> input, fftw_complex* channel_input, uint64_t vector_size) {
    //uint64_t degree = input.size();
    uint64_t degree = 1;
    fftw_complex* ptr_conv_result;
    fftw_complex* ptr_inverse_transform_result;
    std::vector<fftw_complex*> ptr_fftw_out(degree, nullptr);
    fftw_plan plan;
    uint64_t uli, ulj;
    double* ptr_double;
    std::vector<std::vector<std::complex<double> > > complex_input(degree, std::vector<std::complex<double> >(vector_size, (std::complex<double>)0.0));
    std::vector<std::complex<double> > conv_result(vector_size, (std::complex<double>)0.0);
    std::vector<std::complex<double> > complex_channel_input(vector_size, (std::complex<double>)0.0);

/* prepare data area including input data */
    for(uli = 0; uli < degree; uli++) ptr_fftw_out[uli] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    ptr_conv_result = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);

    ptr_inverse_transform_result = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    if(ptr_inverse_transform_result == nullptr) {
        std::cerr << "<variable node update> Can not allocate memory." << std::endl;
        exit(-1);
    }

/* clear output memory area */
    for(uli = 0; uli < degree; uli++) {
        for(ulj = 0; ulj < vector_size; ulj++) {
            ptr_fftw_out[uli][ulj][0] = (double)0.0;
            ptr_fftw_out[uli][ulj][1] = (double)0.0;
        }
    } 

    for(ulj = 0; ulj < vector_size; ulj++) {
        ptr_inverse_transform_result[ulj][0] = (double)0.0;
        ptr_inverse_transform_result[ulj][1] = (double)0.0;
    }

/* execute discrete FFT */
    for(uli = 0; uli < degree; uli++) {
        plan = fftw_plan_dft_1d(vector_size, input[uli], ptr_fftw_out[uli], FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }

/* transform fftw_complex to std::complex */
    for(uli = 0; uli < degree; uli++) {
        for(ulj = 0; ulj < vector_size; ulj++) complex_input[uli][ulj] = density::transform_to_complex(ptr_fftw_out[uli][ulj]);
    }
    for(ulj = 0; ulj < vector_size; ulj++) complex_channel_input[ulj] = density::transform_to_complex(channel_input[ulj]);

/* compute product of fourier transform results */
    for(uli = 0; uli < vector_size; uli++) conv_result[uli] = complex_input[0][uli];
    for(uli = 1; uli < degree; uli++) {
        for(ulj = 0; ulj < vector_size; ulj++) conv_result[ulj] *= complex_input[uli][ulj];
    }
    //std::memcpy(ptr_conv_result, ptr_fftw_out[0], sizeof(fftw_complex)*vector_size);
    //for(uli = 1; uli < degree; uli++) {
    //    for(ulj = 0; ulj < vector_size; ulj++) {
    //        ptr_conv_result[ulj][0] = ptr_conv_result[ulj][0]*ptr_fftw_out[uli][ulj][0] - ptr_conv_result[ulj][1]*ptr_fftw_out[uli][ulj][1];
    //        ptr_conv_result[ulj][1] = ptr_conv_result[ulj][1]*ptr_fftw_out[uli][ulj][0] + ptr_conv_result[ulj][0]*ptr_fftw_out[uli][ulj][1];
    //    }
    //}

    for(ulj = 0; ulj < vector_size; ulj++) conv_result[ulj] *= complex_channel_input[ulj];
    //for(ulj = 0; ulj < vector_size; ulj++) {
    //    ptr_conv_result[ulj][0] = ptr_conv_result[ulj][0]*channel_input[ulj][0] - ptr_conv_result[ulj][1]*channel_input[ulj][1];
    //    ptr_conv_result[ulj][1] = ptr_conv_result[ulj][1]*channel_input[ulj][0] + ptr_conv_result[ulj][0]*channel_input[ulj][1];
    //}
    for(ulj = 0; ulj < vector_size; ulj++) {
        ptr_conv_result[ulj][0] = conv_result[ulj].real();
        ptr_conv_result[ulj][1] = conv_result[ulj].imag();
    }

    plan = fftw_plan_dft_1d(vector_size, ptr_conv_result, ptr_inverse_transform_result, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for(ulj = 0; ulj < vector_size; ulj++) {
        ptr_inverse_transform_result[ulj][0] = ptr_inverse_transform_result[ulj][0]/(double)vector_size;
    }
    ptr_double = (double *)fftw_malloc(sizeof(double)*vector_size);
    if(ptr_double == nullptr) {
        std::cerr << "<variable node update> Can not allocate memory." << std::endl;
        exit(-1);
    }
    /* copy real part */
    for(uli = 0; uli < vector_size; uli++) ptr_double[uli] = ptr_inverse_transform_result[uli][0];
/* free memory */
    for(uli = 0; uli < degree; uli++) {
        fftw_free(ptr_fftw_out[uli]);
    }
    fftw_free(ptr_conv_result);
    fftw_free(ptr_inverse_transform_result);

    return(ptr_double);
}
}