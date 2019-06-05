#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <random>
#include <ctime>
std::default_random_engine generator;
std::uniform_real_distribution<double> distr(0.0, 1.0);
double erand48(unsigned short*) {
	return distr(generator);
}

#include <omp.h>
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif // STB_IMAGE_IMPLEMENTATION

using std::cout;
using std::endl;
using std::min;
using std::max;
const double PI = acos(-1);
const double eps = 1e-6;
const double INF = 1 << 20;

const double min_p[3] = { 100, 100, 100 };
const double max_p[3] = { 10000, 10000, 10000 };

enum Refl_t { DIFF, SPEC, REFR };

int gamma_trans(double x) { return int(.5 + 255 * pow(x < 0 ? 0 : x > 1 ? 1 : x, 1 / 2.2)); }
double sqr(double x) { return x * x; }

