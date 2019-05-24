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
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
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
#include "plot.hpp"

fftw_plan p1,p2,p3;
fftw_complex *_a, *_b, *_c, *_A, *_B, *_C;
double *tmp1, *tmp2, *tmp3, *tmp4;
std::vector<std::vector<int> > R_tbl(vector_size, std::vector<int>(vector_size, 0));
extern degree_distribution load_degree_distribution(std::string);

namespace density {

void clip_pdf(double* x) /* PDF clipping function */
{
  int i;
  double sum;

  sum = 0.0;
  for (i = lower_bound + upper_bound; i < upper_bound + half_lower_bound; i++) {
    sum += x[i];
    x[i] = 0.0;
  }
  x[half_lower_bound + upper_bound] +=sum;
  sum = 0.0;
  for (i = upper_bound + half_upper_bound + 1; i < (int64_t)vector_size; i++) {
    sum += x[i];
    x[i] = 0.0;
  }
  x[half_upper_bound + upper_bound] +=sum;
}

double llr_pdf(double x, double channel_var) /* PDF of LLR (AWGN channel) */
{
  double var = 4.0/channel_var;
  double mean = 2.0/channel_var;
  return (1.0/sqrt(2.0*pi*var)) * exp(-(x-mean)*(x-mean)/(2.0*var));
}

void quantize_pdf(double* a, double channel_var) /* The function makes quatized version of a given PDF function.  */
{
  int i;
  double l,r;
  double vl,vr;

  for (i = half_lower_bound; i < half_upper_bound; i++) {
    l = delta*(double)i - 0.5*delta;
    r = delta*(double)i + 0.5*delta;
    vl = llr_pdf(l, channel_var);
    vr = llr_pdf(r, channel_var);
    a[i + upper_bound] = (vl + vr)*0.5;
  }
  for (i = lower_bound; i < half_lower_bound; i++) a[i + upper_bound] = 0.0;
  for (i = half_upper_bound + 1; i < upper_bound; i++) a[i + upper_bound] = 0.0;
}

double error_prob(double* a)
{
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < upper_bound; i++) {
    sum +=  a[i];
  }
  sum += 0.5*a[upper_bound];
  return(sum*delta);
}

void normalize_pdf(double* r)
{
  int i;
  double sum;

  sum = 0.0;			/* normalization */

  for (i = 0; i < (int64_t)vector_size; i++) sum += r[i]*delta;
  for (i = 0; i < (int64_t)vector_size; i++) r[i]/= sum;
}

void conv(double* a, double* b, double* r) /* convolution of two PDFs */
{
  int i;

  for (i = 0; i < (int64_t)vector_size; i++) {
    _a[i][0] = a[i];
    _a[i][1] = 0.0;
    _b[i][0] = b[i];
    _b[i][1] = 0.0;
  }
  fftw_execute(p1);
  fftw_execute(p2);
  for (i = 0; i < (int64_t)vector_size; i++) {
    _C[i][0] = _A[i][0]*_B[i][0]-_A[i][1]*_B[i][1];
    _C[i][1] = _A[i][0]*_B[i][1]+_A[i][1]*_B[i][0];
  }
  fftw_execute(p3);
  for (i = 0; i < (int64_t)vector_size; i++) {
    r[i] = _c[((vector_size - 1)/2 + i)%vector_size][0]*delta/(double)vector_size;
  }
  clip_pdf(r);
}

void delta_func(double* a) /* delta function */
{
  int i;
  for (i = 0; i < (int64_t)vector_size; i++) {
    a[i] = 0.0;
  }
  a[upper_bound] = 1.0/delta;
}

void naive_convpower(double* a, double* b, int n)
{
  int i;
  for (i = 0; i < (int64_t)vector_size; i++) {
    b[i] = a[i];
  }
  for (i = 1; i < n; i++) {
    conv(a,b,b);
  }
}

void convpower(double* a, double* b, int n)
{
  int i;
  int w;
  int bit;

  if (n < 6) {
    naive_convpower(a, b, n);
    return;
  }
  for (i = 0; i < (int64_t)vector_size; i++) tmp1[i] = a[i];
  delta_func(b);
  w = 0;
  i = 0;
  while(true) {
    bit = (n >> i) & 1;
    if (bit == 1) conv(b, tmp1, b);
    if (bit == 1) w += (1<<i);
    if (w == n) break;
    i++;
    conv(tmp1, tmp1, tmp1);
  }
}

double R(double a, double b)
{
  return 2.0 * atanh(tanh(0.5*a)*tanh(0.5*b));
}

void init(void) /* initialization */
{
  int i,j;
  int idx = 0;
  double a,b,v;

  printf("#vector_size      = %" PRIu64 "\n", vector_size);
  printf("#lower_bound      = %" PRId64 "\n", lower_bound);
  printf("#upper_bound      = %" PRId64 "\n", upper_bound);
  printf("#half_lower_bound = %" PRId64 "\n", half_lower_bound);
  printf("#half_upper_bound = %" PRId64 "\n", half_upper_bound);

  _a = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
  _b = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
  _c = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
  _A = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
  _B = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);
  _C = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*vector_size);

  tmp1 = (double *)fftw_malloc(sizeof(double)*vector_size);
  tmp2 = (double *)fftw_malloc(sizeof(double)*vector_size);
  tmp3 = (double *)fftw_malloc(sizeof(double)*vector_size);
  tmp4 = (double *)fftw_malloc(sizeof(double)*vector_size);
	
  for (i = half_lower_bound; i <= half_upper_bound; i++) {
    for (j = half_lower_bound; j <= half_upper_bound; j++) {
      a = delta*i;
      b = delta*j;
      v = R(a,b);
      if (v >= 0.5*delta) idx = (int)floor(v/delta + 0.5);
      if (v <= -0.5*delta) idx = (int)ceil(v/delta - 0.5);
      if ((v < 0.5*delta) && (v > -0.5*delta)) idx = 0;
      if (v > half_upper_bound*delta) idx = half_upper_bound;
      if (v < half_lower_bound*delta) idx = half_lower_bound;
      R_tbl[i + upper_bound][j + upper_bound] = idx;
    }
  }
  p1 = fftw_plan_dft_1d(vector_size, _a, _A, FFTW_FORWARD, FFTW_ESTIMATE);  
  p2 = fftw_plan_dft_1d(vector_size, _b, _B, FFTW_FORWARD, FFTW_ESTIMATE);  
  p3 = fftw_plan_dft_1d(vector_size, _C, _c, FFTW_BACKWARD, FFTW_ESTIMATE);  
}

void apply_R(double* a, double* b, double* c)
{
  int i,j;

  for (i = 0; i < (int64_t)vector_size; i++) tmp1[i] = 0.0;
  for (i = half_lower_bound + upper_bound; i <= half_upper_bound + upper_bound; i++) {
    for (j = half_lower_bound + upper_bound; j <= half_upper_bound + upper_bound; j++) {
      tmp1[R_tbl[i][j] + upper_bound ] += a[i]*b[j];
    }
  }
  for (i = 0; i < (int64_t)vector_size; i++) c[i] = tmp1[i]*delta;
}

void Rpower(double* a, double* b, int n)
{
  int i;
  for (i = 0; i <(int64_t)vector_size; i++) b[i] = a[i];
  for (i = 1; i < n; i++) apply_R(a,b,b);
}

void add_pdf(double* a, double* b, double weight, uint64_t vector_size) {
  uint64_t uli;

  for(uli = 0; uli < vector_size; uli++) a[uli] += weight*b[uli];
}

void clear_pdf(double* pdf, uint64_t vector_size) {
  for(uint64_t uli = 0; uli < vector_size; uli++) pdf[uli] = (double)0.0;
}

void evolution(uint64_t vector_size, degree_distribution degdist, double channel_var, Plot* plot_app) {
    double* x;
    double probability, probability_old;
    double *P_lambda, *P_c2m, *P_m2c, *f0, *f1;
    double *tmpA;
    double *tmpB;
    double *tmpC;
    uint64_t iteration;
    int64_t i;
    uint64_t uli;

    x = (double *)fftw_malloc(sizeof(double)*vector_size);
    if(x == nullptr) {
        std::cerr << "Can not allocate memory." << std::endl;
        exit(-1);
    }
    for(i = 0; i < (int64_t)vector_size; i++) x[i] = (double)(i - upper_bound)*delta;

    init();

    P_lambda = (double *)fftw_malloc(sizeof(double)*vector_size);
    P_c2m = (double *)fftw_malloc(sizeof(double)*vector_size);
    P_m2c = (double *)fftw_malloc(sizeof(double)*vector_size);
    f0 = (double *)fftw_malloc(sizeof(double)*vector_size);
    f1= (double *)fftw_malloc(sizeof(double)*vector_size);
    tmpA = (double *)fftw_malloc(sizeof(double)*vector_size);
    tmpB = (double *)fftw_malloc(sizeof(double)*vector_size);
    tmpC = (double *)fftw_malloc(sizeof(double)*vector_size);

	/* initialization of LLR density */
    quantize_pdf(P_lambda, channel_var);
    normalize_pdf(P_lambda);

	/* initialization of message(c->m) density */
    delta_func(P_c2m);

    probability_old = 1.0;

    for(iteration = 0; iteration < 1000; iteration++) {
	    /* variable node operation */
      clear_pdf(P_m2c, vector_size);
      clear_pdf(tmpC, vector_size);
      for(uli = 0; uli < degdist.first.size(); uli++) {
        clear_pdf(tmpA, vector_size);
        convpower(P_c2m, tmpA, degdist.first[uli].first);
        add_pdf(tmpC, tmpA, degdist.first[uli].second, vector_size);
      }
      conv(P_lambda, tmpC, P_m2c);
      normalize_pdf(P_m2c);

      /* update graph */
      plot_app->updateCurve(x, P_m2c, vector_size);
      plot_app->emitSignal();
        
		/* BER computation */
      probability = error_prob(P_m2c);
      std::printf("# %" PRIu64 " %16.12f\n", iteration, probability);

      if(fabs(probability-probability_old) < 1e-5) break; /* break condition */

      probability_old = probability;
      clear_pdf(P_c2m, vector_size);
		/* check node operation */ 
      for(uli = 0; uli < degdist.second.size(); uli++) {
        Rpower(P_m2c, tmpC, degdist.second[uli].first);
        normalize_pdf(tmpC);
        add_pdf(P_c2m, tmpC, degdist.second[uli].second, vector_size);
        clear_pdf(tmpC, vector_size);
      }
      normalize_pdf(P_c2m);
    }
    fftw_free(x);
    fftw_free(P_lambda);
    fftw_free(P_c2m);
    fftw_free(P_m2c);
    fftw_free(f0);
    fftw_free(f1);
    fftw_free(tmpA);
    fftw_free(tmpB);
    fftw_free(tmpC);
    fftw_free(_a);
    fftw_free(_b);
    fftw_free(_c);
    fftw_free(_A);
    fftw_free(_B);
    fftw_free(_C);
    fftw_free(tmp1);
    fftw_free(tmp2);
    fftw_free(tmp3);
    fftw_free(tmp4);
}

}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    Plot* probability_plot = new Plot();
    double* x;
    std::string curve_name;
    int64_t i;
    double sigma;
    degree_distribution degdist;

    if(argc != 3) {
        std::cerr << "Usage: ./density <degree distribution file name> <sigma>" << std::endl;
        return(EXIT_FAILURE);
    }
    sigma = std::stod(argv[2]);
    degdist = load_degree_distribution(std::string(argv[1]));

    x = (double *)fftw_malloc(sizeof(double)*vector_size);
    if(x == nullptr) {
        std::cerr << "Can not allocate memory." << std::endl;
        return(EXIT_FAILURE);
    }
    for(i = 0; i < (int64_t)vector_size; i++) x[i] = (double)(i - upper_bound)*delta;

    // draw probability density function 
    probability_plot->resize(800, 600);
    probability_plot->setTitle("probability density function");
    curve_name = R"(probability function)";
    probability_plot->setupCurve(curve_name);
    probability_plot->setAxisScale(2, x[0], x[vector_size - 1], 0);
    probability_plot->updateAxes();
    QObject::connect(probability_plot, SIGNAL(emitSignal()), probability_plot, SLOT(replot()));
    probability_plot->show();

     std::thread evo_thread(density::evolution, vector_size, degdist, sigma*sigma, probability_plot);
     evo_thread.detach();

    return(app.exec());
}