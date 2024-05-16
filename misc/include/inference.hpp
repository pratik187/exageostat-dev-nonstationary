#ifndef INFERENCE_HPP
#define INFERENCE_HPP
#include <bits/stdc++.h>
using namespace std;

extern "C" int predict(double *x, double *y, double* Z_values, double startx, double endx, double starty, double endy, double startz , double endz , bool z_user_range , int N);

#endif
