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

extern "C" {
    #include <fftw3.h>
}

namespace density {

fftw_complex* update_variable_node(std::vector<fftw_complex*> input, fftw_complex* channel_input, uint64_t vector_size) {
    uint64_t degree = input.size();
    fftw_complex* ptr_conv_result;
    fftw_complex* ptr_inverse_transform_result;
    std::vector<fftw_complex*> ptr_fftw_out(degree, nullptr);
    fftw_plan plan;
    uint64_t uli, ulj;

/* prepare data area including input data */
    for(uli = 0; uli < degree; uli++) ptr_fftw_out[uli] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    ptr_conv_result = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    ptr_inverse_transform_result = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);

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
    }

/* compute product of fourier transform results */
    for(uli = 0; uli < degree; uli++) {
        for(ulj = 0; ulj < vector_size; ulj++) {
            ptr_conv_result[ulj][0] = ptr_conv_result[ulj][0]*ptr_fftw_out[uli][ulj][0] - ptr_conv_result[ulj][1]*ptr_fftw_out[uli][ulj][1];
            ptr_conv_result[ulj][1] = ptr_conv_result[ulj][1]*ptr_fftw_out[uli][ulj][0] + ptr_conv_result[ulj][0]*ptr_fftw_out[uli][ulj][1];
        }
    }

    for(ulj = 0; ulj < vector_size; ulj++) {
        ptr_conv_result[ulj][0] = ptr_conv_result[ulj][0]*channel_input[ulj][0] - ptr_conv_result[ulj][1]*channel_input[ulj][1];
        ptr_conv_result[ulj][1] = ptr_conv_result[ulj][1]*channel_input[ulj][0] + ptr_conv_result[ulj][0]*channel_input[ulj][1];
    }

    plan = fftw_plan_dft_1d(vector_size, ptr_conv_result, ptr_inverse_transform_result, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

/* free memory */
    fftw_destroy_plan(plan);
    for(uli = 0; uli < degree; uli++) fftw_free(ptr_fftw_out[uli]);
    fftw_free(ptr_conv_result);
    return(ptr_inverse_transform_result);
}
}