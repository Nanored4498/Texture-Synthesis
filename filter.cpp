#include <math.h>
#include "filter.hpp"

#define uchar unsigned char
#define Complex std::complex<double>
#define VC std::vector<Complex>
#define VVC std::vector<VC>

void average(uchar* in, uchar* out, int size, int W, int H) {
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			int miX = std::max(0, i - size), maX = std::min(W-1, i + size);
			int miY = std::max(0, j - size), maY = std::min(H-1, j + size);
			int d = (maX - miX + 1) * (maY - miY + 1);
			for(int c = 0; c < 3; c++) {
				double sum = 0;
				for(int x = miX; x <= maX; x ++)
					for(int y = miY; y <= maY; y++)
						sum += in[(x + W*y) * 3 + c];
				out[(i + W*j) * 3 + c] = std::max(0.0, std::min(255.0, sum / d));
			}
		}
	}
}

void fft(VC& x) {
    int N = x.size();
    if(N <= 1) return;
	VC even, odd;
	for(int i = 0; i < N; i++) {
		if(i % 2 == 0) even.push_back(x[i]);
		else odd.push_back(x[i]);
	}
	fft(even);
	fft(odd);
	for(int k = 0; k < N/2; k++) {
    	Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

void ifft(VC& x) {
	int N = x.size();
	for(int i = 0; i < N; i++) x[i] = std::conj(x[i]);
	fft(x);
	for(int i = 0; i < N; i++) x[i] = std::conj(x[i]) / ((double) N);
}

void pad(VC& x) {
	int p = 1, N = x.size();
	while(p < N) p *= 2;
	while(N++ < p) x.push_back(0); 
}

void get_cols(uchar* in, VVC &cols, int c, int W, int H) {
	for(int i = 0; i < W; i++) {
		VC co;
		for(int j = 0; j < H; j++) co.push_back(in[(i + j*W) * 3 + c]);
		cols.push_back(co);
	}
}

void cols_to_rows(VVC &cols, VVC &rows, int W, int H) {
	rows.clear();
	for(int i = 0; i < H; i++) {
		VC r;
		for(int j = 0; j < W; j++) r.push_back(cols[j][i]);
		rows.push_back(r);
	}
}

void apply_filter_1D(VC &data, VC& filter) {
	int N = data.size();
	pad(data);
	int N2 = data.size();
	fft(data);
	for(int i = 0; i < N2; i++)
		data[i] *= filter[i];
	ifft(data);
	data.resize(N);
}

void gauss_ifilter(int N, VC& filter, double sigma) {
	int p = 1;
	while(p < N) p *= 2;
	N = p;
	double s = 0;
	for(int i = 0; i < N; i++) {
		double v = i < N/2 ? i : N - i;
		v = std::exp(-v*v / sigma);
		s += v;
		filter.push_back(v);
	}
	for(int i = 0; i < N; i++) filter[i] /= s;
}

void gauss_filter(int N, VC &filter, double sigma) {
	gauss_ifilter(N, filter, sigma);
	fft(filter);
}

void dgauss_filter(int N, VC &filter, double sigma) {
	gauss_ifilter(N, filter, sigma);
	N = filter.size();
	double s = 0;
	for(int i = 0; i < N; i++) {
		filter[i] *= - (i < N / 2 ? i : i - N) / sigma;
		s += std::abs(filter[i]);
	}
	for(int i = 0; i < N; i++) filter[i] *= 2 / s;
	fft(filter);
}

void rows_to_uchar(uchar* out, VVC rows, int W, int H, int c) {
	for(int i = 0; i < W; i ++)
		for(int j = 0; j < H; j++)
			out[(i + j*W)*3 + c] = std::max(0.0, std::min(255.0, rows[j][i].real()));
}

void double_filter_1D(uchar* in, uchar* out, int W, int H, VC &filter_col, VC &filter_row) {
	for(int i = 0; i < 3*W*H; i++) out[i] = in[i];
	for(int c = 0; c < 3; c++) {
		VVC cols, rows;
		get_cols(in, cols, c, W, H);
		for(int i = 0; i < W; i++) apply_filter_1D(cols[i], filter_col);
		cols_to_rows(cols, rows, W, H);
		for(int i = 0; i < H; i++) apply_filter_1D(rows[i], filter_row);
		rows_to_uchar(out, rows, W, H, c);
	}
}

void gauss(uchar* in, uchar* out, double sigma, int W, int H) {
	VC filter_col, filter_row;
	gauss_filter(H, filter_col, sigma*sigma);
	gauss_filter(W, filter_row, sigma*sigma);
	double_filter_1D(in, out, W, H, filter_col, filter_row);
}

void dgauss(uchar* in, uchar* out, double sigma, int W, int H) {
	VC filter_col, filter_row;
	dgauss_filter(H, filter_col, sigma*sigma);
	dgauss_filter(W, filter_row, sigma*sigma);
	double_filter_1D(in, out, W, H, filter_col, filter_row);
}

void bilateral(uchar* in, uchar* out, int size, double sigma1, double sigma2, int W, int H) {
	sigma1 *= sigma1;
	sigma2 *= sigma2;
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			int miX = std::max(0, i - size), maX = std::min(W-1, i + size);
			int miY = std::max(0, j - size), maY = std::min(H-1, j + size);
			for(int c = 0; c < 3; c++) {
				double colij = in[(i + W*j) * 3 + c];
				double sum = 0, w = 0;
				for(int x = miX; x <= maX; x ++) {
					for(int y = miY; y <= maY; y++) {
						double col = in[(x + W*y) * 3 + c];
						double dx = x - i, dy = y - j;
						double a = std::exp( - abs(colij - col) / sigma1 - (dx*dx + dy*dy) / sigma2);
						sum += col * a;
						w += a;
					}
				}
				out[(i + W*j) * 3 + c] = std::max(0.0, std::min(255.0, sum / w));
			}
		}
	}
}