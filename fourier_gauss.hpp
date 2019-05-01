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
#include <cmath>
#include <vector>

namespace density {
    double get_llr_pdf(double x, double mean, double sigma) {
        double rtnv;
        double variance = sigma*sigma;

        rtnv = exp(-0.5*(x - mean)*(x - mean)/variance)/sqrt(2.0*M_PI*variance);
        return(rtnv);
    }

    std::vector<double> quantize_pdf(double sigma) {
        double left_variable, right_variable;
        double left_pdf_value, right_pdf_value;
        std::vector<double> rtnv(vector_size, 0.0);
        int64_t i;
        double mean = ((double)half_upper_bound)*delta;

        for(i = half_lower_bound; i < upper_bound + 1; i++) {
            left_variable = delta*i - 0.5*delta;
            right_variable = delta*i + 0.5*delta;
            left_pdf_value = get_llr_pdf(left_variable, mean, sigma);
            right_pdf_value = get_llr_pdf(right_variable, mean, sigma);
            rtnv[i + upper_bound] = 0.5*(left_pdf_value + right_pdf_value);
        }
        //for(i = lower_bound; i < half_lower_bound; i++) rtnv[i + upper_bound] = 0.0;
        //for(i = half_upper_bound + 1; i < upper_bound; i++) rtnv[i] = 0.0;
        return(rtnv);
    }

    void normalize_pdf(double* pdf, uint64_t vector_size) {
        std::vector<double> rtnv(vector_size, 0.0);
        uint64_t uli;
        double sum = 0;

        for(uli = 0; uli < vector_size; uli++) sum += pdf[uli];

        for(uli = 0; uli < vector_size; uli++) pdf[uli] = pdf[uli]/sum;
    }

    double get_error_probability(double* pdf) {
        double sum = 0.0;
        uint64_t uli;

        for(uli = 0; uli < upper_bound; uli++) sum += pdf[uli];
        sum += 0.5*pdf[upper_bound];
        std::cout << "sum = " << sum << std::endl;
        sum = sum*delta;
        return(sum);
    }

    std::vector<double> clip_pdf(std::vector<double> pdf) {
        std::vector<double> rtnv(vector_size, 0.0);
        double sum = 0.0;
        int64_t i;
		for(i = lower_bound + upper_bound; i < upper_bound + half_lower_bound; i++) {
			sum += pdf[i];
			rtnv[i] = 0.0;
		}
		rtnv[half_lower_bound + upper_bound] +=sum;
		sum = 0.0;
		for(i = upper_bound + half_upper_bound + 1; i < (int64_t)vector_size; i++) {
			sum += pdf[i];
			rtnv[i] = 0.0;
		}
		rtnv[half_upper_bound + upper_bound] += sum;
        return(rtnv);
	}

    std::vector<double> get_delta_function(void) {
        std::vector<double> rtnv(vector_size, 0.0);
        uint64_t uli;

        for(uli = 0; uli < vector_size; uli++) rtnv[uli] = 0.0;
        rtnv[upper_bound + half_upper_bound] = 1.0/delta;
        return(rtnv);
    }

    double get_absolute_value(double a, double b) {
        double absv = sqrt(a*a + b*b);
        return(absv);
    }

    int64_t quasi_gallager(int64_t x, int64_t y, uint64_t vector_size) {
        double a, b;
        double z = 0.0;
       uint64_t rtnv;

        a = 0.5*((double)x)*delta;
        b = 0.5*((double)y)*delta;

        z = 2.0*atanh(tanh(a)*tanh(b));
        rtnv = (int64_t)std::round(z/delta);
        return(rtnv);
    }

    double* get_two_way_tanh_map(double* a, double* b, uint64_t vector_size) {
        int64_t i, j, k;
        double* ptr_output = (double *)fftw_malloc(sizeof(double)*vector_size);

        for(i = 0; i < (int64_t)vector_size; i++) {
            ptr_output[i] = 0.0;
        }

        for(i = lower_bound; i < upper_bound + 1; i++) {
            for(j = lower_bound; j < upper_bound + 1; j++) {
                k = quasi_gallager(i, j, vector_size);
                ptr_output[k + upper_bound] += a[i + upper_bound]*b[j + upper_bound];
            }
        }
        return(ptr_output);
    }

    double* convolute_pdf(double* a, double* b, uint64_t vector_size) {
        fftw_complex* fftw_complex_a;
        fftw_complex* fftw_complex_b;
        fftw_complex* fftw_complex_c;
        uint64_t uli;
        double* rtnv;
        uint64_t complex_vector_size = 2*(vector_size/2) + 1;
        fftw_plan plan;

        fftw_complex_a = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
        plan = fftw_plan_dft_r2c_1d(vector_size, a, fftw_complex_a, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_complex_b = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
        plan = fftw_plan_dft_r2c_1d(vector_size, b, fftw_complex_b, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_complex_c = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*complex_vector_size);
        for(uli = 0; uli < complex_vector_size; uli++) {
            fftw_complex_c[uli][0] = fftw_complex_a[uli][0]*fftw_complex_b[uli][0] - fftw_complex_a[uli][1]*fftw_complex_b[uli][1];
            fftw_complex_c[uli][1] = fftw_complex_a[uli][1]*fftw_complex_b[uli][0] + fftw_complex_a[uli][0]*fftw_complex_b[uli][1];
        }
        rtnv = (double *)fftw_malloc(sizeof(double)*vector_size);
        plan = fftw_plan_dft_c2r_1d(complex_vector_size, fftw_complex_c, rtnv, FFTW_ESTIMATE);
        fftw_execute(plan);
        for(uli = 0; uli < vector_size; uli++) rtnv[uli] = rtnv[uli]/vector_size;
        fftw_destroy_plan(plan);
        fftw_free(fftw_complex_a);
        fftw_free(fftw_complex_b);
        fftw_free(fftw_complex_c);
        return(rtnv);
    }
}