// -*- c++ -*-

// 1. implement "fixed" generator
// 2. implement discrete generator
// 3. implement combine generator? 

#ifndef GENERATOR_H
#define GENERATOR_H

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#include "config.h"

#include <string>
#include <vector>
#include <utility>

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <random>
#include <cmath>

#include "log.h"
#include "util.h"
#include "include/nlohmann/json.hpp"

using json = nlohmann::json;

// Generator syntax:
//
// \d+ == fixed
// n[ormal]:mean,sd
// e[xponential]:lambda
// p[areto]:scale,shape
// g[ev]:loc,scale,shape
// fb_value, fb_key, fb_rate

class Generator {
public:
  Generator() {}
  //  Generator(const Generator &g) = delete;
  //  virtual Generator& operator=(const Generator &g) = delete;
  virtual ~Generator() {}

  virtual double generate(double U = -1.0) = 0;
  virtual void set_lambda(double lambda) {DIE("set_lambda() not implemented");}
protected:
  std::string type;
};

class Fixed : public Generator {
public:
  Fixed(double _value = 1.0) : value(_value) { D("Fixed(%f)", value); }
  virtual double generate(double U = -1.0) { return value; }
  virtual void set_lambda(double lambda) {
    if (lambda > 0.0) value = 1.0 / lambda;
    else value = 0.0;
  }

private:
  double value;
};

class Uniform : public Generator {
public:
  Uniform(double _scale) : scale(_scale) { D("Uniform(%f)", scale); }

  virtual double generate(double U = -1.0) {
    if (U < 0.0) U = drand48();
    return scale * U;
  }

  virtual void set_lambda(double lambda) {
    if (lambda > 0.0) scale = 2.0 / lambda;
    else scale = 0.0;
  }

private:
  double scale;
};

class Normal : public Generator {
public:
  Normal(double _mean = 1.0, double _sd = 1.0) : mean(_mean), sd(_sd) {
    D("Normal(mean=%f, sd=%f)", mean, sd);
  }

  virtual double generate(double U = -1.0) {
    if (U < 0.0) U = drand48();
    double V = U; // drand48();
    double N = sqrt(-2 * log(U)) * cos(2 * M_PI * V);
    return mean + sd * N;
  }

  virtual void set_lambda(double lambda) {
    if (lambda > 0.0) mean = 1.0 / lambda;
    else mean = 0.0;
  }

private:
  double mean, sd;
};

class Exponential : public Generator {
public:
  Exponential(double _lambda = 1.0) : lambda(_lambda) {
    D("Exponential(lambda=%f)", lambda);
  }

  virtual double generate(double U = -1.0) {
    if (lambda <= 0.0) return 0.0;
    if (U < 0.0) U = drand48();
    return -log(U) / lambda;
  }

  virtual void set_lambda(double lambda) { this->lambda = lambda; }

private:
  double lambda;
};

class NH_Exponential : public Generator {
  public:
  NH_Exponential(double _lambda = 1.0) : lambda(_lambda) {
      D("NH_Exponential(lambda=%f)", lambda);
      load_A_histogram("./asset/amplitude_counts.json", histogram);
      load_phi_histogram("./period_counts.json", histogram);
      load_w_histogram("./phase_counts.json", histogram);
    }
  
    // 速率函数 lambda(t) = A * sin(w * t + phi) + b
    // 累积强度函数 Lambda(t) = integral(lambda(u), 0, t)，这个case下是解析解
    double cumulative_intensity(double t, double A, double w, double phi, double b) {
      return (A / w) * (-std::cos(w * t + phi) + std::cos(phi)) + b * t;
    }

    void load_A_histogram(const std::string& filename) {
      std::ifstream file(filename);
      if (!file.is_open()) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return;
      }

      json j;
      file >> j;
      file.close();

      double total_count = 0.0;
      for (auto& [bin, count] : j.items()) {
        A_histogram[bin] = count;
          total_count += A_histogram[bin];
      }

      // Normalize counts to probabilities
      for (auto& [bin, count] : A_histogram) {
          double count_double = count;
          A_histogram[bin] = count_double / total_count;
      }
    }

    void load_phi_histogram(const std::string& filename) {
      std::ifstream file(filename);
      if (!file.is_open()) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return;
      }

      json j;
      file >> j;
      file.close();

      double total_count = 0.0;
      for (auto& [bin, count] : j.items()) {
        A_histogram[bin] = count;
          total_count += phi_histogram[bin];
      }

      // Normalize counts to probabilities
      for (auto& [bin, count] : phi_histogram) {
          double count_double = count;
          phi_histogram[bin] = count_double / total_count;
      }
    }

    void load_w_histogram(const std::string& filename) {
      std::ifstream file(filename);
      if (!file.is_open()) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return;
      }

      json j;
      file >> j;
      file.close();

      double total_count = 0.0;
      for (auto& [bin, count] : j.items()) {
        A_histogram[bin] = count;
          total_count += w_histogram[bin];
      }

      // Normalize counts to probabilities
      for (auto& [bin, count] : w_histogram) {
          double count_double = count;
          w_histogram[bin] = count_double / total_count;
      }
    }

    // 逆累积强度函数：求解 Lambda(T) = x 得到 T
    double inverse_cumulative_intensity(double x, double A, double   w, double phi, double b) {
      double t = 0.0; // 初始猜测
      double epsilon = 1e-6; // 精度
      double step = 0.01; // 步长

      // 使用数值方法求解 Lambda(t) = x
      while (cumulative_intensity(t, A, w, phi, b) < x) {
          t += step;
      }
      return t;
    }

    // 生成非齐次指数分布的随机数
    double generate_nonhomogeneous_exponential(double A, double w, double phi, double b) {
      double u = drand48();
      double x = -std::log(1.0 - u); // 转换为指数分布的随机数
      return inverse_cumulative_intensity(x, A, w, phi, b); // 通过逆累积强度函数得到 T
    }

    virtual double generate(double U = -1.0) {
      if (lambda <= 0.0) return 0.0;
      adjust_A();
      adjust_phi();
      adjust_w();
      return generate_nonhomogeneous_exponential(A, w, phi, b);
    }
  
    virtual void set_lambda(double lambda) { this->lambda = lambda; }
  
  private:
    double lambda;
    double A = 2;
    double w = 0.5;
    double phi = 0;
    double b = 3;

    std::map<std::string, double> A_histogram;
    std::map<std::string, double> w_histogram;
    std::map<std::string, double> phi_histogram;

    void adjust_A() {
      static std::random_device rd;
      static std::mt19937 gen(rd());
      static std::uniform_real_distribution<> dis(0.0, 1.0);

      double random_value = dis(gen);
      double cumulative_probability = 0.0;

      for (const auto& [bin, prob] : A_histogram) {
          cumulative_probability += prob;
          if (random_value <= cumulative_probability) {
              // Parse the bin range and set A to the midpoint of the bin
              size_t dash_pos = bin.find('-');
              double lower_bound = std::stod(bin.substr(0, dash_pos));
              double upper_bound = std::stod(bin.substr(dash_pos + 1));
              A = (lower_bound + upper_bound) / 2.0;
              break;
          }
      }
    }

    void adjust_w() {
      static std::random_device rd;
      static std::mt19937 gen(rd());
      static std::uniform_real_distribution<> dis(0.0, 1.0);

      double random_value = dis(gen);
      double cumulative_probability = 0.0;

      for (const auto& [bin, prob] : A_histogram) {
          cumulative_probability += prob;
          if (random_value <= cumulative_probability) {
              // Parse the bin range and set A to the midpoint of the bin
              size_t dash_pos = bin.find('-');
              double lower_bound = std::stod(bin.substr(0, dash_pos));
              double upper_bound = std::stod(bin.substr(dash_pos + 1));
              w = (lower_bound + upper_bound) / 2.0;
              break;
          }
      }
    }

    void adjust_phi() {
      static std::random_device rd;
      static std::mt19937 gen(rd());
      static std::uniform_real_distribution<> dis(0.0, 1.0);

      double random_value = dis(gen);
      double cumulative_probability = 0.0;

      for (const auto& [bin, prob] : A_histogram) {
          cumulative_probability += prob;
          if (random_value <= cumulative_probability) {
              // Parse the bin range and set A to the midpoint of the bin
              size_t dash_pos = bin.find('-');
              double lower_bound = std::stod(bin.substr(0, dash_pos));
              double upper_bound = std::stod(bin.substr(dash_pos + 1));
              phi = (lower_bound + upper_bound) / 2.0;
              break;
          }
      }
    }
  };

class GPareto : public Generator {
public:
  GPareto(double _loc = 0.0, double _scale = 1.0, double _shape = 1.0) :
    loc(_loc), scale(_scale), shape(_shape) {
    assert(shape != 0.0);
    D("GPareto(loc=%f, scale=%f, shape=%f)", loc, scale, shape);
  }

  virtual double generate(double U = -1.0) {
    if (U < 0.0) U = drand48();
    return loc + scale * (pow(U, -shape) - 1) / shape;
  }

  virtual void set_lambda(double lambda) {
    if (lambda <= 0.0) scale = 0.0;
    else scale = (1 - shape) / lambda - (1 - shape) * loc;
  }

private:
  double loc /* mu */;
  double scale /* sigma */, shape /* k */;
};

class GEV : public Generator {
public:
  GEV(double _loc = 0.0, double _scale = 1.0, double _shape = 1.0) :
    e(1.0), loc(_loc), scale(_scale), shape(_shape) {
    assert(shape != 0.0);
    D("GEV(loc=%f, scale=%f, shape=%f)", loc, scale, shape);
  }

  virtual double generate(double U = -1.0) {
    return loc + scale * (pow(e.generate(U), -shape) - 1) / shape;
  }

private:
  Exponential e;
  double loc /* mu */, scale /* sigma */, shape /* k */;
};

class Discrete : public Generator {
public:
  ~Discrete() { delete def; }
  Discrete(Generator* _def = NULL) : def(_def) {
    if (def == NULL) def = new Fixed(0.0);
  }

  virtual double generate(double U = -1.0) {
    double Uc = U;
    if (pv.size() > 0 && U < 0.0) U = drand48();

    double sum = 0;
 
    for (auto p: pv) {
      sum += p.first;
      if (U < sum) return p.second;
    }

    return def->generate(Uc);
  }

  void add(double p, double v) {
    pv.push_back(std::pair<double,double>(p, v));
  }

private:
  Generator *def;
  std::vector< std::pair<double,double> > pv;
};

class KeyGenerator {
public:
  KeyGenerator(Generator* _g, double _max = 10000) : g(_g), max(_max) {}
  std::string generate(uint64_t ind) {
    uint64_t h = fnv_64(ind);
    double U = (double) h / ULLONG_MAX;
    double G = g->generate(U);
    int keylen = MAX(round(G), floor(log10(max)) + 1);
    char key[256];
    snprintf(key, 256, "%0*" PRIu64, keylen, ind);

    //    D("%d = %s", ind, key);
    return std::string(key);
  }
private:
  Generator* g;
  double max;
};

Generator* createGenerator(std::string str);
Generator* createFacebookKey();
Generator* createFacebookValue();
Generator* createFacebookIA();

#endif // GENERATOR_H
