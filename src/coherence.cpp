#include "coherence.hpp"
#include <iostream>
#include <cmath>
#include "stb_image.h"
#include "stb_image_write.h"

void appendCoherence(VVP &C, uchar* E, int m, int m2, int h) {
	std::vector<std::vector<VI>> Ns(m);
	int kn = std::min(7, m / h / 2);
	if(kn % 2 == 0) kn ++;
	for(int x = 0; x < m; x++) {
		Ns[x] = std::vector<VI>(m);
		for(int y = 0; y < m; y++)
			getNeighb(x, y, kn, E, m2, Ns[x][y], h);
	}
	double md = pow(0.35*m/kn, 2.0);
	int i = 0;
	std::cerr << "Coherence (" << m << ") [";
	for(int sp = 0; sp < (m>>2); sp++) std::cerr << " ";
	std::cerr << "]\rCoherence (" << m << ") [";
	for(int y = 0; y < m; y++) {
		#pragma omp parallel for
		for(int x = 0; x < m; x++) {
			VI Na = Ns[x][y];
			Pix best;
			int dist = 1e9;
			i = x+m*y;
			for(int x2 = 0; x2 < m; x2++) {
				for(int y2 = 0; y2 < m; y2++) {
					bool cont = false;
					for(Pix &p : C[i]) {
						int dx = abs(p.x - x2), dy = abs(p.y - y2);
						if(pow(std::min(dx, m2-dx), 2) + pow(std::min(dy, m2-dy), 2) < md) {
							cont = true;
							break;
						}
					}
					if(cont) continue;
					int d = distNeighb(Na, Ns[x2][y2], kn);
					if(d < dist) {
						dist = d;
						best.x = x2;
						best.y = y2;
					}
				}
			}
			C[i].push_back(best);
		}
		if(y % 4 == 3) std::cerr << "=" << std::flush;
	}
	std::cerr << "]" << std::endl;
}

void initCoherence(VVP &C, int m) {
	for(int y = 0; y < m; y++)
		for(int x = 0; x < m; x++)
			C.push_back({{x, y}});
}

void loadCoherence(VVP &C, const char* filename) {
	int m, c;
	uchar* im = stbi_load(filename, &m, &m, &c, 3);
	int i = 0;
	for(int y = 0; y < m; y++)
		for(int x = 0; x < m; x++)
			C[i++].push_back({im[3*(x+y*m)], im[3*(x+y*m)+1]});
	free(im);
}

void writeCoherence(VVP &C, int k, const char* filename, int m) {
	uchar* coh = new uchar[3*m*m]();
	for(int i = 0; i < m*m; i++)
		coh[3*i] = C[i][k].x, coh[3*i+1] = C[i][k].y;
	stbi_write_png(filename, m, m, 3, coh, 0);
	delete[] coh;
}