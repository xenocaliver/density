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
#include "common.hpp"
#include "transform_complex.hpp"

namespace density {
void clip_pdf(double* pdf, uint64_t vector_size) {
	int64_t i;
	double sum = 0.0;

	for (i = lower_bound + upper_bound; i < upper_bound + half_lower_bound; i++) {
		sum += pdf[i];
		pdf[i] = 0.0;
	}
	pdf[half_lower_bound + upper_bound] +=sum;
	sum = 0.0;
	for (i = upper_bound + half_upper_bound + 1; i < (int64_t)vector_size; i++) {
		sum += pdf[i];
	    pdf[i] = 0.0;
	}
	pdf[half_upper_bound + upper_bound] += sum;
}
double* update_variable_node(std::vector<double*> input, fftw_complex* channel_input, uint64_t vector_size) {
    uint64_t degree = input.size();
    uint64_t complex_vector_size = vector_size/2 + 1;
    fftw_complex* ptr_conv_result;
    double* ptr_inverse_transform_result;
    std::vector<fftw_complex*> ptr_fftw_out(degree, nullptr);
    fftw_plan plan;
    uint64_t uli, ulj;
    std::vector<std::vector<std::complex<double> > > complex_input(degree, std::vector<std::complex<double> >(complex_vector_size, (std::complex<double>)0.0));
    std::vector<std::complex<double> > conv_result(complex_vector_size, (std::complex<double>)0.0);
    std::vector<std::complex<double> > complex_channel_input(complex_vector_size, (std::complex<double>)0.0);

/* prepare data area including input data */
    for(uli = 0; uli < degree; uli++)  {
        ptr_fftw_out[uli] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*complex_vector_size);
        if(ptr_fftw_out[uli] == nullptr) {
            std::cerr << "Can not allocate memory." << std::endl;
            exit(-1);
        }
    }

    ptr_conv_result = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*complex_vector_size);
    if(ptr_conv_result == nullptr) {
        std::cerr << "Can not allocate memory." << std::endl;
        exit(-1);
    }

    ptr_inverse_transform_result = (double *)fftw_malloc(sizeof(double)*vector_size);
    if(ptr_inverse_transform_result == nullptr) {
        std::cerr << "<variable node update> Can not allocate memory." << std::endl;
        exit(-1);
    }

/* clear output memory area */
    for(uli = 0; uli < degree; uli++) {
        for(ulj = 0; ulj < complex_vector_size; ulj++) {
            ptr_fftw_out[uli][ulj][0] = (double)0.0;
            ptr_fftw_out[uli][ulj][1] = (double)0.0;
        }
    } 

    for(ulj = 0; ulj < vector_size; ulj++) {
        ptr_inverse_transform_result[ulj] = (double)0.0;
    }

/* execute discrete FFT */
    for(uli = 0; uli < degree; uli++) {
        plan = fftw_plan_dft_r2c_1d(vector_size, input[uli], ptr_fftw_out[uli], FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }

/* transform fftw_complex to std::complex */
    for(uli = 0; uli < degree; uli++) {
        for(ulj = 0; ulj < complex_vector_size; ulj++) complex_input[uli][ulj] = density::transform_to_complex(ptr_fftw_out[uli][ulj]);
    }
    for(ulj = 0; ulj < complex_vector_size; ulj++) complex_channel_input[ulj] = density::transform_to_complex(channel_input[ulj]);

/* compute product of fourier transform results */
    for(uli = 0; uli < complex_vector_size; uli++) conv_result[uli] = complex_input[0][uli];
    for(uli = 1; uli < degree; uli++) {
        for(ulj = 0; ulj < complex_vector_size; ulj++) conv_result[ulj] *= complex_input[uli][ulj];
    }

    for(ulj = 0; ulj < complex_vector_size; ulj++) conv_result[ulj] *= complex_channel_input[ulj];

    ptr_conv_result = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*complex_vector_size);
    if(ptr_conv_result == nullptr) {
        std::cerr << "Can not allocate memory." << std::endl;
        exit(-1);
    }
    for(ulj = 0; ulj < complex_vector_size; ulj++) {
        ptr_conv_result[ulj][0] = conv_result[ulj].real();
        ptr_conv_result[ulj][1] = conv_result[ulj].imag();
    }

    plan = fftw_plan_dft_c2r_1d(vector_size, ptr_conv_result, ptr_inverse_transform_result, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for(ulj = 0; ulj < vector_size; ulj++) {
        ptr_inverse_transform_result[ulj] = ptr_inverse_transform_result[ulj]*delta/(double)vector_size;
    }

/* free memory */
    for(uli = 0; uli < degree; uli++) {
        fftw_free(ptr_fftw_out[uli]);
    }
    fftw_free(ptr_conv_result);

    return(ptr_inverse_transform_result);
}

}