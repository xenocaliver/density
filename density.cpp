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
#include <string>
#include <boost/property_tree/json_parser.hpp>

#include <QtCore>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>

extern "C" {
#include <fftw3.h>
}
#include "common.hpp"
#include "fourier_gauss.hpp"
#include "plot.hpp"

namespace density {

extern double* update_variable_node(std::vector<double*>, fftw_complex*, uint64_t);
extern double* update_check_node(std::vector<double*>, uint64_t);

void evolution(std::vector<double*> variable_node_inputs, fftw_complex* channel_dft_input, std::vector<double*> check_node_inputs, uint64_t vector_size, Plot* plot_app) {
    uint64_t iteration;
    double* ptr_double;
    double* ptr_result;
    double probability;
    double* x;
    uint64_t uli, ulj;
    int64_t i;
    std::vector<double*> v(variable_node_inputs.size() + 1, nullptr);

    x = (double *)fftw_malloc(sizeof(double)*vector_size);
    if(x == nullptr) {
        std::cerr << "Can not allocate memory." << std::endl;
        exit(-1);
    }
    for(i = 0; i < (int64_t)vector_size; i++) x[i] = ((double)(i - upper_bound))*delta;

    for(uli = 0; uli < v.size(); uli++) {
        v[uli] = (double *)fftw_malloc(sizeof(double)*vector_size);
        if(v[uli] == nullptr) {
            std::cerr << "Can not allocate memory." << std::endl;
            exit(-1);
        }
    }

    for(iteration = 0; iteration < 40; iteration++) {
    /*variable node calculation */
        ptr_double = density::update_variable_node(variable_node_inputs, channel_dft_input, vector_size);
    
    /* copy probability distribution */
        for(uli = 0; uli < 5; uli++) std::memcpy(check_node_inputs[uli], ptr_double, sizeof(double)*vector_size);

    /* free memrmoy */
        fftw_free(ptr_double);

    /* check node calculation */
        ptr_double = density::update_check_node(check_node_inputs, vector_size);
        density::normalize_pdf(ptr_double, vector_size);

    /* copy memory */
        for(uli = 0; uli < 2; uli++) {
            for(ulj = 0; ulj < vector_size; ulj++) {
                variable_node_inputs[uli][ulj] = ptr_double[ulj];
            }
        }
   /* get joint probability density */
        for(uli = 0; uli < v.size(); uli++) std::memcpy(v[uli], ptr_double, sizeof(double)*vector_size);
        fftw_free(ptr_double);
        ptr_result = density::update_variable_node(v, channel_dft_input, vector_size);
        normalize_pdf(ptr_result, vector_size);

    /* update graph */
        plot_app->updateCurve(x, ptr_result, vector_size);
        plot_app->emitSignal();

    /* caculating joint probability */
        probability = density::get_error_probability(ptr_result);
        std::cout << std::setw(8) << probability << std::endl;
        fftw_free(ptr_result);
    }
    for(uli = 0; uli < 2; uli++) fftw_free(variable_node_inputs[uli]);
    for(uli = 0; uli < 5; uli++) fftw_free(check_node_inputs[uli]);
    for(uli = 0; uli < v.size(); uli++) fftw_free(v[uli]);
}

}

int main(int argc, char* argv[]) {
    uint64_t uli;
    std::vector<double*> variable_node_inputs(2, nullptr);
    std::vector<double*> check_node_inputs(5, nullptr);
    double* channel_input = (double *)fftw_malloc(sizeof(double)*vector_size);
    fftw_complex* channel_dft_input = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*complex_vector_size);
    std::vector<double> initial_distribution;
    std::vector<double> delta_distribution;
    double sigma = 0.8;
    double* ptr_input;
    fftw_plan plan;
    double* ptr_double;
    QApplication app(argc, argv);
    Plot* delta_plot = new Plot();
    Plot* awgn_plot = new Plot();
    Plot* probability_plot = new Plot();
    double* x;
    std::string curve_name;
    int64_t i;

    x = (double *)fftw_malloc(sizeof(double)*vector_size);
    if(x == nullptr) {
        std::cerr << "Can not allocate memory." << std::endl;
        return(EXIT_FAILURE);
    }
    for(i = 0; i < (int64_t)vector_size; i++) x[i] = ((double)(i - upper_bound))*delta;

    for(uli = 0; uli < 2; uli++) {
        ptr_input = (double *)fftw_malloc(sizeof(double)*vector_size);
        if(ptr_input == nullptr) {
            std::cerr << "Can not allocate memory." << std::endl;
            exit(-1);
        }
        variable_node_inputs[uli] = ptr_input;
    }
    for(uli = 0; uli < 5; uli++) {
        ptr_double = (double *)fftw_malloc(sizeof(double)*vector_size);
        if(ptr_double == nullptr) {
            std::cerr << "Can not allocate memory." << std::endl;
            exit(-1);
        }
        check_node_inputs[uli] = ptr_double;
    }

    /* initial distribution */
    initial_distribution = density::quantize_pdf(sigma);
    density::normalize_pdf(initial_distribution.data(), initial_distribution.size());
    for(uli = 0; uli < vector_size; uli++) {
        channel_input[uli] = initial_distribution[uli];
    }

    plan = fftw_plan_dft_r2c_1d(vector_size, channel_input, channel_dft_input, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    delta_distribution = density::get_delta_function();
    density::normalize_pdf(delta_distribution.data(), delta_distribution.size());
    for(uli = 0; uli < vector_size; uli++) {
        variable_node_inputs[0][uli] = delta_distribution[uli];
    }

    // draw input delta function
    delta_plot->resize(800, 600);
    delta_plot->setTitle("delta distribution");
    curve_name = R"(delta function)";
    delta_plot->setupCurve(curve_name);
    delta_plot->setAxisScale(0, 0.0, 1.2, 0);
    delta_plot->setAxisScale(2, x[0], x[vector_size - 1] + 1.0, 0);
    delta_plot->updateAxes();
    delta_plot->plotCurve(x, delta_distribution.data(), vector_size);
    delta_plot->show();

    // draw input AWGN pdf
    awgn_plot->resize(800, 600);
    awgn_plot->setTitle("AWGN noise distribution");
    curve_name = R"(AWGN)";
    awgn_plot->setupCurve(curve_name);
    awgn_plot->setAxisScale(2, x[0], x[vector_size - 1] + 1.0, 0);
    awgn_plot->updateAxes();
    awgn_plot->plotCurve(x, initial_distribution.data(), vector_size);
    awgn_plot->show();

    // draw probability density function 
    probability_plot->resize(800, 600);
    probability_plot->setTitle("probability density function");
    curve_name = R"(probability function)";
    probability_plot->setupCurve(curve_name);
    probability_plot->setAxisScale(2, x[0], x[vector_size - 1] + 1.0, 0);
    probability_plot->updateAxes();
    QObject::connect(probability_plot, SIGNAL(emitSignal()), probability_plot, SLOT(replot()));
    probability_plot->show();

    std::memcpy(variable_node_inputs[1], variable_node_inputs[0], sizeof(double)*vector_size);

     std::thread evo_thread(density::evolution, variable_node_inputs, channel_dft_input, check_node_inputs, vector_size, probability_plot);
     evo_thread.detach();

    return(app.exec());
}