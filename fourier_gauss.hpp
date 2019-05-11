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

    void clip_pdf(double* pdf) {
        double sum = 0.0;
        int64_t i;
		for(i = lower_bound + upper_bound; i < upper_bound + half_lower_bound; i++) {
			sum += pdf[i];
			pdf[i] = 0.0;
		}
		pdf[half_lower_bound + upper_bound] +=sum;
		sum = 0.0;
		for(i = upper_bound + half_upper_bound + 1; i < (int64_t)vector_size; i++) {
			sum += pdf[i];
			pdf[i] = 0.0;
		}
		pdf[half_upper_bound + upper_bound] += sum;
	}

    std::vector<double> get_delta_function(void) {
        std::vector<double> rtnv(vector_size, 0.0);
        uint64_t uli;

        for(uli = 0; uli < vector_size; uli++) rtnv[uli] = 0.0;
        rtnv[upper_bound + half_upper_bound] = 1.0/delta;
        return(rtnv);
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

    void update_tanh_map(int64_t i, int64_t j, double* a, double* b, double* output, uint64_t vector_size) {
        int64_t k;

        k = quasi_gallager(i, j, vector_size);
        output[k + upper_bound] += a[i + upper_bound]*b[j + upper_bound];
    }

    void wrapper_update_tanh_map(std::vector<int64_t> bound, double* a, double* b, double* output, uint64_t vector_size) {
        int64_t i, j;

        for(i = bound[0]; i < bound[1] + 1; i++) {
            for(j = bound[2]; j < bound[3] + 1; j++) {
                update_tanh_map(i, j, a, b, output, vector_size);
            }
        }
    }

    double* get_two_way_tanh_map(double* a, double* b, uint64_t vector_size) {
        int64_t i, j;
        std::vector<double*> ptr_output(4, nullptr);
        std::vector<std::vector<int64_t> > bound(4, std::vector<int64_t>(4, 0));

        for(i = 0; i < 4; i++) {
            ptr_output[i] = (double *)fftw_malloc(sizeof(double)*vector_size);
            if(ptr_output[i] == nullptr) {
                std::cerr << "Can not allocate memory." << std::endl;
                exit(-1);
            }
            for(j = 0; j < (int64_t)vector_size; j++) ptr_output[i][j] = (double)0.0;
        }

        bound[0][0] = lower_bound;
        bound[0][1] = 0;
        bound[0][2] = lower_bound;
        bound[0][3] = 0; 

        bound[1][0] = 0;
        bound[1][1] = upper_bound;
        bound[1][2] = lower_bound;
        bound[1][3] = 0; 

        bound[2][0] = lower_bound;
        bound[2][1] = 0;
        bound[2][2] = 0;
        bound[2][3] = upper_bound; 

        bound[3][0] = 0;
        bound[3][1] = upper_bound;
        bound[3][2] = 0;
        bound[3][3] = upper_bound; 

        try {
            std::thread t1(wrapper_update_tanh_map, bound[0], a, b, ptr_output[0], vector_size);
            std::thread t2(wrapper_update_tanh_map, bound[1], a, b, ptr_output[1], vector_size);
            std::thread t3(wrapper_update_tanh_map, bound[2], a, b, ptr_output[2], vector_size);
            std::thread t4(wrapper_update_tanh_map, bound[3], a, b, ptr_output[3], vector_size);
            t4.join();
            t3.join();
            t2.join();
            t1.join();
        } catch(std::exception &e) {
            std::cerr << e.what() << std::endl;
            exit(-1);
        }
        for(uint64_t uli = 0; uli < vector_size; uli++) ptr_output[0][uli] += ptr_output[1][uli] + ptr_output[2][uli] + ptr_output[3][uli];
        for(uint64_t uli = 1; uli < 4; uli++) fftw_free(ptr_output[uli]);
        return(ptr_output[0]);
    }
}