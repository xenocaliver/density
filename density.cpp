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
#include <cstring>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>

extern "C" {
#include <fftw3.h>
}
#include "common.hpp"
#include "fourier_gauss.hpp"
#include "QWindows.hpp"

namespace density {

extern double* update_variable_node(std::vector<fftw_complex*>, fftw_complex*, uint64_t);
extern double* update_check_node(std::vector<double*>, uint64_t);

double integrate_real_part(fftw_complex* pdf_real, uint64_t vector_size) {
    uint64_t uli;
    double rtnv = 0.0;

    for(uli = 0; uli < vector_size; uli++) {
        rtnv += pdf_real[uli][0];
    }
    return(rtnv);
}
}

int main(int argc, char* argv[]) {
    uint64_t uli;
    std::vector<fftw_complex*> variable_node_inputs(2, nullptr);
    std::vector<double*> check_node_inputs(5, nullptr);
    fftw_complex* channel_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    fftw_complex* channel_dft_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
    double probability;
    uint64_t iteration;
    std::vector<double> initial_distribution;
    std::vector<double> delta_distribution;
    double sigma = 0.5;
    fftw_complex* ptr_input;
    fftw_plan plan;
    double* ptr_double;
    double* ptr_tmp;
    std::ofstream ofs;
    QApplication app(argc, argv);

    ViewPDF *v = new ViewPDF();
    v->show();
    app.exec();

    for(uli = 0; uli < 2; uli++) {
        ptr_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
        variable_node_inputs[uli] = ptr_input;
    }
    for(uli = 0; uli < 5; uli++) {
        ptr_double = (double *)fftw_malloc(sizeof(double)*vector_size);
        check_node_inputs[uli] = ptr_double;
    }

    /* initial distribution */
    initial_distribution = density::quantize_pdf(sigma);
    density::normalize_pdf(initial_distribution.data(), initial_distribution.size());
    ofs.open("initial.txt");
    if(!ofs) {
        std::cerr << "Can not open file: initial.txt" << std::endl; 
        return(EXIT_FAILURE);
    }
    for(uli = 0; uli < vector_size; uli++) {
        channel_input[uli][0] = initial_distribution[uli];
        channel_input[uli][1] = (double)0.0;
        ofs << channel_input[uli][0] << std::endl;
    }
    ofs.close();

    plan = fftw_plan_dft_1d(vector_size, channel_input, channel_dft_input, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    delta_distribution = density::get_delta_function();
    density::normalize_pdf(delta_distribution.data(), delta_distribution.size());
    ofs.open("delta.txt");
    if(!ofs) {
        std::cerr << "Can not open file: delta.txt" << std::endl;
        return(EXIT_FAILURE);
    }
    for(uli = 0; uli < vector_size; uli++) {
        variable_node_inputs[0][uli][0] = delta_distribution[uli];
        variable_node_inputs[0][uli][1] = (double)0.0;
        ofs << variable_node_inputs[0][uli][0] << std::endl;
    }
    ofs.close();

    std::memcpy(variable_node_inputs[1], variable_node_inputs[0], sizeof(fftw_complex)*vector_size);

    for(iteration = 0; iteration < 5; iteration++) {
    /*variable node calculation */
        ptr_double = density::update_variable_node(variable_node_inputs, channel_dft_input, vector_size);
        density::normalize_pdf(ptr_double, vector_size);

    /* caculating joint probability */
        probability = density::get_error_probability(ptr_double);
        std::cout << std::setw(8) << probability << std::endl;
        ofs.open("test" + std::to_string(iteration) + ".txt");
        for(uli = 0; uli < vector_size; uli++) ofs << ptr_double[uli] << std::endl;
        ofs.close();
        fftw_free(ptr_double);
        return(EXIT_SUCCESS);
    /* copy probability distribution */
        for(uli = 0; uli < 5; uli++) std::memcpy(check_node_inputs[uli], ptr_double, sizeof(double)*vector_size);

    /* free memrmoy */
        fftw_free(ptr_double);

    /* check node calculation */
        ptr_double = density::update_check_node(check_node_inputs, vector_size);
        density::normalize_pdf(ptr_double, vector_size);

        std::memcpy(variable_node_inputs[0], ptr_double, sizeof(double)*vector_size);
        std::memcpy(variable_node_inputs[1], ptr_double, sizeof(double)*vector_size);
        fftw_free(ptr_double);
    }

    for(uli = 0; uli < 2; uli++) fftw_free(variable_node_inputs[uli]);
    for(uli = 0; uli < 5; uli++) fftw_free(check_node_inputs[uli]);
    return(EXIT_SUCCESS);
}