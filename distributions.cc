#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <random>
#include <cmath>

#include "distributions.h"
#include "log.h"

const char* distributions[] =
  { "uniform", "exponential", "zipfian", "latest", NULL };

distribution_t get_distribution(const char *name) {
  for (int i = 0; distributions[i] != NULL; i++)
    if (!strcmp(distributions[i], name))
      return (distribution_t) i;
  return (distribution_t) -1;
}

double generate_normal(double mean, double sd) {
  double U = drand48();
  double V = drand48();
  double N = sqrt(-2 * log(U)) * cos(2 * M_PI * V);
  return mean + sd * N;
}

double generate_poisson(double lambda) {
  if (lambda <= 0.0) return 0;
  double U = drand48();
  return -log(U)/lambda;
}

double generate_uniform(double lambda) {
  if (lambda <= 0.0) return 0;
  return 1.0 / lambda;
}

// // 速率函数 lambda(t) = A * sin(w * t + phi) + b
// // 累积强度函数 Lambda(t) = integral(lambda(u), 0, t)，这个case下是解析解
// double cumulative_intensity(double t, double A, double w, double phi, double b) {
//     return (A / w) * (-std::cos(w * t + phi) + std::cos(phi)) + b * t;
// }

// // 逆累积强度函数：求解 Lambda(T) = x 得到 T
// double inverse_cumulative_intensity(double x, double A, double w, double phi, double b) {
//     double t = 0.0; // 初始猜测
//     double epsilon = 1e-6; // 精度
//     double step = 0.01; // 步长

//     // 使用数值方法求解 Lambda(t) = x
//     while (cumulative_intensity(t, A, w, phi, b) < x) {
//         t += step;
//     }
//     return t;
// }

// // 生成非齐次指数分布的随机数
// double generate_nonhomogeneous_exponential(double A, double w, double phi, double b) {
//     double u = drand48();
//     double x = -std::log(1.0 - u); // 转换为指数分布的随机数
//     return inverse_cumulative_intensity(x, A, w, phi, b); // 通过逆累积强度函数得到 T
// }