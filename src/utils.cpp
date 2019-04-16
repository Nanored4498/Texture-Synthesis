#include "utils.hpp"
#include "stb_image_write.h"
#include <cmath>

void getNeighb_S(int x, int y, int k, Pix* S, int W, int H, uchar* E, int m, VI &N) {
	N.clear();
	int dk = (k-1) / 2;
	x += W;
	y += H;
	for(int i = -dk; i <= dk; i++) {
		for(int j = -dk; j <= dk; j++) {
			int indS = ((x + i) % W) + ((y + j) % H) * W;
			Pix u = S[indS];
			int indE = 3 * (u.x + u.y*m);
			N.push_back(E[indE]);
			N.push_back(E[indE+1]);
			N.push_back(E[indE+2]);
		}
	}
}

void getNeighb(int x, int y, int k, uchar* E, int m, VI &N, int h) {
	N.clear();
	int dk = (k-1) / 2;
	for(int i = -dk; i <= dk; i++) {
		for(int j = -dk; j <= dk; j++) {
			int indE = 3 * (((m + x + h*i) % m) + ((m + y + h*j) % m) * m);
			N.push_back(E[indE]);
			N.push_back(E[indE+1]);
			N.push_back(E[indE+2]);
		}
	}
}

bool is_tor(uchar* im, int m) {
	int d = 0;
	for(int c = 0; c < 3; c++) {
		for(int i = 0; i < m; i++) {
			int diffx = im[3*i*m+c] - im[3*(m-1+m*i)+c];
			int diffy = im[3*i+c] - im[3*(i+m*(m-1))+c];
			d += diffx*diffx + diffy*diffy;
		}
	}
	return d < 1444 * 2 * m * 3;
}

uchar* torrify(uchar* im, int m) {
	int m2 = m << 1;
	uchar* res = new uchar[4*m2*m2];
	for(int c = 0; c < 3; c++) {
		for(int x = 0; x < m; x++) {
			for(int y = 0; y < m; y++) {
				uchar v = im[3*(x+m*y)+c];
				int x2 = m2-1-x, y2 = m2-1-y;
				res[3*(x+m2*y)+c] = v;
				res[3*(x+m2*y2)+c] = v;
				res[3*(x2+m2*y)+c] = v;
				res[3*(x2+m2*y2)+c] = v;
			}
		}
	}
	return res;
}

uchar* square(uchar* im, int W, int H) {
	int m = std::min(W, H);
	uchar* res = new uchar[3*m*m];
	for(int c = 0; c < 3; c++)
		for(int x = 0; x < m; x++)
			for(int y = 0; y < m; y++)
				res[3*(x+m*y)+c] = im[3*(x+W*y)+c];
	return res;
}

void save(Pix* S, int W, int H, uchar* E, int m, const char* name) {
	int size = W*H;
	uchar* res = new uchar[3*size];
	for(int i = 0; i < size; i++) {
		int u = 3*(S[i].x + m*S[i].y);
		res[3*i] = E[u];
		res[3*i+1] = E[u+1];
		res[3*i+2] = E[u+2];
	}
	stbi_write_png(name, W, H, 3, res, 0);
	delete[] res;
}

void saveS(Pix* S, int W, int H, int m, const char* name) {
	int size = W*H;
	double scale = 256.0 / m;
	uchar* res = new uchar[3*size]();
	for(int i = 0; i < size; i++)
		res[3*i] = scale*S[i].x, res[3*i+1] = scale*S[i].y;
	stbi_write_png(name, W, H, 3, res, 0);
	delete[] res;
}

const static double MEAN_FILTER[3][3] = {{0.05, 0.15, 0.05}, {0.15, 0.2, 0.15}, {0.05, 0.15, 0.05}};
void save_smooth(Pix* S, int W, int H, uchar* E, int m, const char* name) {
	int size = W*H;
	uchar* res = new uchar[3*size];
	for(int i = 0; i < size; i++) {
		int x = i % W, y = i / W;
		bool mean = (x > 0 && x < W-1 && y > 0 && y < H-1)
				&& (abs(S[i].x - S[i-1].x) + abs(S[i].y - S[i-1].y) > 4
					|| abs(S[i].x - S[i+1].x) + abs(S[i].y - S[i+1].y) > 4
					|| abs(S[i].x - S[i-W].x) + abs(S[i].y - S[i-W].y) > 4
					|| abs(S[i].x - S[i+W].x) + abs(S[i].y - S[i+W].y) > 4);
		if(mean) {
			for(int c = 0; c < 3; c++) {
				double sum = 0;
				for(int dx = -1; dx <= 1; dx ++) {
					for(int dy = -1; dy <= 1; dy ++) {
						int j = i + dx + W*dy;
						int u = 3*(S[j].x + m*S[j].y);
						sum += MEAN_FILTER[dx+1][dy+1] * E[u+c];
					}
				}
				res[3*i+c] = sum;
			}
		} else {
			int u = 3*(S[i].x + m*S[i].y);
			res[3*i] = E[u];
			res[3*i+1] = E[u+1];
			res[3*i+2] = E[u+2];
		}
	}
	stbi_write_png(name, W, H, 3, res, 0);
	delete[] res;
}