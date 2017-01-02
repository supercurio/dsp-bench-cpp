/*
 *  Copyright Francois Simond 2017
 */

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>

using std::cout;
using std::vector;

struct Biquad {
  double b0;
  double b1;
  double b2;
  double a1;
  double a2;

  double x1;
  double x2;
  double y1;
  double y2;
};

Biquad peakEqCoefs(double Fs, double f0, double Q, double dBgain) {
  double A = pow(10.0, dBgain / 40.0);
  double omega = 2 * M_PIl * f0 / Fs;

  double alpha = sin(omega) / (2 * Q);

  Biquad bq;

  double a0 = 1 + alpha / A;

  bq.b0 = (1 + alpha * A) / a0;
  bq.b1 = (-2 * cos(omega)) / a0;
  bq.b2 = (1 - alpha * A) / a0;
  bq.a1 = bq.b1;
  bq.a2 = (1 - alpha / A) / a0;

  return bq;
}

void printBiquad(Biquad bq) {
  cout << "Biquad coefficients:\n";

  cout << std::setprecision(17);

  cout << "b0=" << bq.b0 << ",\n";
  cout << "b1=" << bq.b1 << ",\n";
  cout << "b2=" << bq.b2 << ",\n";
  cout << "a1=" << bq.a1 << ",\n";
  cout << "a2=" << bq.a2 << ",\n";
}

void whiteNoise(vector<double> &buffer) {
  unsigned int seed = time(0);
  for (unsigned int i = 0; i < buffer.size(); i++) {
    buffer[i] = static_cast<double>(rand_r(&seed)) / RAND_MAX;
  }
}

void iir(const vector<double> &in, vector<double> &out, Biquad bq) {
  for (unsigned int i = 0; i < in.size(); i++) {
    out[i] = (bq.b0 * in[i]) + (bq.b1 * bq.x1) + (bq.b2 * bq.x2) -
             (bq.a1 * bq.y1) - (bq.a2 * bq.y2);

    bq.x2 = bq.x1;
    bq.x1 = in[i];

    bq.y2 = bq.y1;
    bq.y1 = out[i];
  }
}

int main(int argc, char **argv) {
  cout << "DSP Bench C++\n";

  size_t buffer_len = 4096;
  if (argc == 2) buffer_len = atoi(argv[1]);

  const int benchLoops = 200000;

  Biquad bq = peakEqCoefs(48000, 200, 2, 6);
  printBiquad(bq);

  vector<double> input;
  input.resize(buffer_len);
  whiteNoise(input);

  vector<double> output;
  output.resize(input.size());

  auto start = std::chrono::steady_clock::now();  // first time_point
  for (int i = 0; i < benchLoops; i++) {
    iir(input, output, bq);
  }
  auto end = std::chrono::steady_clock::now();  // second time_point

  auto nanos =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

  cout << nanos.count() / benchLoops << " ns per loop" << '\n';

  return 0;
}
