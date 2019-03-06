#include <iostream>
#include <vector>
#include <queue>
#include "filter.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace std;

typedef unsigned char uchar;
typedef vector<double> VD;
typedef vector<int> VI;

struct Pix {
	int x, y;
	Pix(int a=0, int b=0): x(a), y(b) {}
};
Pix operator*(int a, const Pix &p) { return {a*p.x, a*p.y}; }
Pix operator+(const Pix &a, const Pix &b) { return {a.x+b.x, a.y+b.y}; }
Pix operator%(const Pix &p, int d) { return {p.x % d, p.y % d}; }
typedef vector<Pix> VP;
typedef vector<VP> VVP;

Pix* upsample(Pix* S, int W, int h, int m) {
	Pix* res = (Pix*) malloc(W*W*4*sizeof(Pix));
	int W2 = W*2;
	h = floor(h*0.5);
	Pix ze = {m-h, m-h}, ri = {h, m-h}, up = {m-h, h}, di = {h, h};
	// Pix ri = {h, 0}, up = {0, h}, di = {h, h};
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < W; j++) {
			// Pix u = 2 * S[i+W*j];
			Pix u = S[i+W*j];
			int p = 2*(i + j*W2);
			res[p] = (u + ze) % m;
			res[p+1] = (u + ri) % m;
			res[p+W2] = (u + up) % m;
			res[p+W2+1] = (u + di) % m;
		}
	}
	return res;
}

unsigned int int_hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}
Pix hp(int p) {
	int a = abs((int) int_hash(p));
	int b = abs((int) int_hash(a+p));
	return {2*(a % 2) - 1, 2*(b % 2) - 1};
}

void jitter(Pix* S, int W, int h, int m, double r) {
	double hr = h*r;
	int W2 = W*W;
	for(int i = 0; i < W2; i++) {
		Pix q = hp(i);
		S[i] = (S[i] + Pix(floor(hr*q.x + 0.5), floor(hr*q.y + 0.5))) % m;
	}
}

void getNeighb_S(int x, int y, int k, Pix* S, int W, uchar* E, int m, VI &N) {
	N.clear();
	int dk = (k-1) / 2;
	x += W;
	y += W;
	for(int i = -dk; i <= dk; i++) {
		for(int j = -dk; j <= dk; j++) {
			int indS = ((x + i) % W) + ((y + j) % W) * W;
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

int distNeighb(VI &a, VI &b, int k) {
	int k2 = k*k*3;
	int res = 0;
	for(int i = 0; i < k2; i++)
		res += pow(a[i] - b[i], 2);
	return res;
}

void correct(Pix* S, int W, int h, int m, VVP &C, uchar* E, double kappa) {
	VI NSp, NEu;
	double phi2 = pow(1 + kappa, 2.0);
	VP sps = {{0, 0}, {1, 1}, {0, 1}, {1, 0}, {1, 1}, {0, 0}, {0, 1}, {1, 0}};
	for(Pix &sp : sps) {
		#pragma omp parallel for private(NSp, NEu)
		for(int px = sp.x; px < W; px += 2) {
			for(int py = sp.y; py < W; py += 2) {
				int dis = 1e9;
				Pix umin;
				getNeighb_S(px, py, 5, S, W, E, m, NSp);
				for(int i = -1; i <= 1; i++) {
					for(int j = -1; j <= 1; j++) {
						int dpx = (W+px+i) % W, dpy = (W+py+j) % W;
						Pix du = (S[dpx+W*dpy] + Pix(m-h*i, m-h*j)) % m;
						int du_ind = du.x + m*du.y;
						for(int k = 0; k < C[du_ind].size(); k++) {
							Pix cu = C[du_ind][k];
							getNeighb(cu.x, cu.y, 5, E, m, NEu, h);
							int d = distNeighb(NSp, NEu, 5);
							if(k > 0) d *= phi2;
							if(d < dis) {
								dis = d;
								umin = cu;
							}
						}
					}
				}
				S[px + W*py] = umin;
			}
		}
	}
}

void save(Pix* S, int ss, uchar* E, int W, int H, const char* name) {
	int s2 = ss*ss;
	uchar* res = (uchar*) malloc(s2*3);
	for(int i = 0; i < s2; i++) {
		int u = 3*(S[i].x + W*S[i].y);
		res[3*i] = E[u];
		res[3*i+1] = E[u+1];
		res[3*i+2] = E[u+2];
	}
	stbi_write_png(name, ss, ss, 3, res, 0);
	delete res;
}
void saveS(Pix* S, int W, const char* name) {
	int W2 = W*W;
	uchar* res = (uchar*) calloc(W2*3, 1);
	for(int i = 0; i < W2; i++)
		res[3*i] = S[i].x, res[3*i+1] = S[i].y;
	stbi_write_png(name, W, W, 3, res, 0);
	delete res;
}

void appendCoherence(VVP &C, uchar* E, int m) {
	vector<vector<VI>> Ns(m);
	for(int x = 0; x < m; x++) {
		Ns[x] = vector<VI>(m);
		for(int y = 0; y < m; y++)
			getNeighb(x, y, 7, E, m, Ns[x][y], 1);
	}
	double md = pow(0.05*m, 2.0);
	int i = 0;
	for(int y = 0; y < m; y++) {
		for(int x = 0; x < m; x++) {
			VI Na = Ns[x][y];
			Pix best;
			int dist = 1e9;
			#pragma omp parallel for
			for(int x2 = 0; x2 < m; x2++) {
				for(int y2 = 0; y2 < m; y2++) {
					bool cont = false;
					for(Pix &p : C[i]) {
						int dx = abs(p.x - x2), dy = abs(p.y - y2);
						if(pow(min(dx, m-dx), 2) + pow(min(dy, m-dy), 2) < md) {
							cont = true;
							break;
						}
					}
					if(cont) continue;
					int d = distNeighb(Na, Ns[x2][y2], 7);
					#pragma omp critical
					if(d < dist) {
						dist = d;
						best = {x2, y2};
					}
				}
			}
			C[i++].push_back(best);
		}
		cerr << "Coherence line: " << y << endl;
	}
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
	delete im;
}

void writeCoherence(VVP &C, int k, const char* filename, int m) {
	uchar* coh = (uchar*) calloc(m*m*3, 1);
	for(int i = 0; i < m*m; i++)
		coh[3*i] = C[i][k].x, coh[3*i+1] = C[i][k].y;
	stbi_write_png(filename, m, m, 3, coh, 0);
	delete coh;
}

Pix* synthesize(int m, VD &r, int c, double kappa, uchar* E) {
	Pix* S = (Pix*) malloc(4*sizeof(Pix));
	uchar* El = (uchar*) malloc(3*m*m);
	int L = log2(m);
	S[0] = {0, 0};
	S[1] = {0, 0};
	S[2] = {0, 0};
	S[3] = {0, 0};
	int W = 2;
	char name[100];
	for(int l = 1; l <= L; l++) {
		int h = 1 << (L-l);
		gauss(E, El, h, m, m);
		Pix* nS = upsample(S, W, h, m);
		delete S;
		S = nS;
		W *= 2;	
		if(l < r.size() && r[l] > 0)
			jitter(S, W, h, m, r[l]);
		if(l > 2) {
			VVP C;
			initCoherence(C, m);
			sprintf(name, "coherence_1_%d.png", l);
			// loadCoherence(C, name);		/* if coherence has been pre_computed then load it */
			// appendCoherence(C, El, m);		/* To compute coherence. May take a while */
			// writeCoherence(C, 1, name, m);		/* To save coherence in a file */
			for(int i = 0; i < c; i++)
				correct(S, W, h, m, C, El, kappa);
		}
		sprintf(name, "out_%d.png", l);
		save(S, W, El, m, m, name);
		sprintf(name, "map_%d.png", l);
		saveS(S, W, name);
		sprintf(name, "pyramid_%d.png", l);
		stbi_write_png(name, m, m, 3, El, 0);
	}
	delete El;
	return S;
}

int main() {
	int W, H, col;
	uchar* E = stbi_load("ims/1.png", &W, &H, &col, 3);
	if (!E) {
		std::cout<<"loading failed"<<std::endl;
		return 1;
	}
	VD r = {0.2, 0.3, 0.6, 0.3, 0.4, 0.1, 0.1, 0.3};
	Pix* S = synthesize(W, r, 2, 0.2, E);
	delete E, S;
	return 0;
}