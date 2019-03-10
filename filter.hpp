#ifndef DEF_FRACTION_H
#define DEF_FRACTION_H

#include <vector>
#include <complex>

#define uchar unsigned char
#define Complex std::complex<double>
#define VC std::vector<Complex>
#define VVC std::vector<VC>

void average(uchar* in, uchar* out, int size, int W, int H);

void fft(VC& x);

void ifft(VC& x);

void gauss(uchar* in, uchar* out, double sigma, int W, int H);

void dgauss(uchar* in, uchar* out, double sigma, int W, int H);

void bilateral(uchar* in, uchar* out, int size, double sigma1, double sigma2, int W, int H);

#undef uchar
#undef Complex
#undef VC
#undef VCC

#endif