#include <iostream>
#include <vector>
#include <queue>
#include "filter.hpp"
#include "stb_image.h"
#include "stb_image_write.h"

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
Pix operator-(const Pix &a, const Pix &b) { return {a.x-b.x, a.y-b.y}; }
Pix operator%(const Pix &p, int d) { return {p.x % d, p.y % d}; }
Pix operator/(const Pix &p, int d) { return {p.x / d, p.y / d}; }
typedef vector<Pix> VP;
typedef vector<VP> VVP;

Pix* upsample(Pix* S, int W, int H, int h, int m) {
	Pix* res = new Pix[W*H*4];
	int W2 = W*2;
	h = floor(h*0.5);
	Pix ze = {m-h, m-h}, ri = {h, m-h}, up = {m-h, h}, di = {h, h};
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
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
	unsigned int a = int_hash(p);
	unsigned int b = int_hash(a+p);
	return {(int) ((b & 1) << 1) - 1, (int) ((b & 1) << 1) - 1};
}

void jitter(Pix* S, int W, int H, int h, int m, double r) {
	double hr = h*r;
	int size = W*H;
	for(int i = 0; i < size; i++) {
		Pix q = hp(i);
		S[i] = (S[i] + Pix(m + floor(hr*q.x + 0.5), m + floor(hr*q.y + 0.5))) % m;
	}
}

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

int distNeighb(VI &a, VI &b, int k) {
	int k2 = k*k*3;
	int res = 0;
	for(int i = 0; i < k2; i++)
		res += pow(a[i] - b[i], 2);
	return res;
}

void correct(Pix* S, int W, int H, int h, int m, int m2, VVP &C, uchar* E, double kappa) {
	VI NSp, NEu;
	double phi2 = pow(1 + kappa, 2.0);
	VP sps = {{0, 0}, {1, 1}, {0, 1}, {1, 0}, {1, 1}, {0, 0}, {0, 1}, {1, 0}};
	for(Pix &sp : sps) {
		#pragma omp parallel for private(NSp, NEu)
		for(int px = sp.x; px < W; px += 2) {
			for(int py = sp.y; py < H; py += 2) {
				int dis = 1e9;
				Pix umin;
				getNeighb_S(px, py, 5, S, W, H, E, m2, NSp);
				for(int i = -1; i <= 1; i++) {
					for(int j = -1; j <= 1; j++) {
						int dpx = (W+px+i) % W, dpy = (H+py+j) % H;
						Pix du = (S[dpx+W*dpy] + Pix(m-h*i, m-h*j)) % m;
						int du_ind = du.x + m*du.y;
						int c_size = C[du_ind].size();
						for(int k = 0; k < c_size; k++) {
							Pix cu = C[du_ind][k];
							getNeighb(cu.x, cu.y, 5, E, m2, NEu, h);
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

void appendCoherence(VVP &C, uchar* E, int m, int h) {
	vector<vector<VI>> Ns(m);
	for(int x = 0; x < m; x++) {
		Ns[x] = vector<VI>(m);
		for(int y = 0; y < m; y++)
			getNeighb(x, y, min(7, m / h), E, m, Ns[x][y], h);
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
	uchar* coh = new uchar[3*m*m]();
	for(int i = 0; i < m*m; i++)
		coh[3*i] = C[i][k].x, coh[3*i+1] = C[i][k].y;
	stbi_write_png(filename, m, m, 3, coh, 0);
	delete[] coh;
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
	int m = min(W, H);
	uchar* res = new uchar[3*m*m];
	for(int c = 0; c < 3; c++)
		for(int x = 0; x < m; x++)
			for(int y = 0; y < m; y++)
				res[3*(x+m*y)+c] = im[3*(x+W*y)+c];
	return res;
}

uchar* downsample(uchar* im, int m, int p) {
	uchar* g = new uchar[3*m*m];
	int h = 1 << p;
	gauss(im, g, h, m, m);
	int m2 = m >> p;
	uchar* res = new uchar[3*m2*m2];
	for(int c = 0; c < 3; c++)
		for(int x = 0; x < m2; x++)
			for(int y = 0; y < m2; y++)
				res[3*(x+m2*y)+c] = g[3*h*(x+m*y)+c];
	delete[] g;
	return res;
}

uchar* magnify(int ml, uchar* Eh, int mh, Pix* S, int W, int H) {
	int mp = pow(2.0, ceil(log2(ml)));
	int Wl = mp*W, Hl = mp*H;
	int Wh = (mh*Wl) / ml, Hh = (mh*Hl) / ml;
	uchar* res = new uchar[3*Wh*Hh];
	double colors[2][2] = {{0, 0}, {0, 0}};
	Pix u;
	for(int c = 0; c < 3; c++) {
		for(int x = 0; x < Wh; x++) {
			for(int y = 0; y < Hh; y++) {
				int xl = (x*ml) / mh, yl = (y*ml) / mh;
				Pix dp = Pix(x*ml - xl*mh, y*ml - yl*mh);
				double fx = (double) dp.x / (double) mh;
				double fy = (double) dp.y / (double) mh;
				dp = dp / ml;
				for(int i = 0; i < 2; i++) {
					for(int j = 0; j < 2; j++) {
						if(xl+i >= Wl || yl+j >= Hl) colors[i][j] = colors[0][0];
						else {
							u = (mh * (S[xl+i + Wl*(yl+j)] - Pix(i, j))) / ml + dp;
							colors[i][j] = Eh[3 * (max(0, u.x) + mh*max(0, u.y)) + c];
						}
					}
				}
				res[3*(x + Wh*y)+c] = (colors[0][0]*(1-fy)+colors[0][1]*fy) * (1-fx)
									+ (colors[1][0]*(1-fy)+colors[1][1]*fy) * fx;
			}
		}
	}
	return res;
}

Pix* synthesize(uchar* E, int m, VD &r, int c, double kappa, int W, int H, bool tor=false) {
	Pix* S = new Pix[W*H];
	uchar* El = new uchar[3*m*m];
	uchar *tE = NULL, *tEl = NULL;
	if(tor) {
		tE = torrify(E, m);
		tEl = new uchar[3*m*m*4];
	}
	int L = ceil(log2(m));
	int r_size = r.size();
	for(int i = 0; i < W*H; i++)
		S[i] = {0, 0};
	if(r_size > 0 && r[0] > 0)
		jitter(S, W, H, m, m, r[0]);
	char name[100];
	for(int l = 1; l <= L; l++) {
		int h = 1 << (L-l);
		if(h == 1) {
			delete[] El;
			El = E;
		} else gauss(E, El, h, m, m);
		if(tor) {
			if(h == 1) {
				delete tEl;
				tEl = tE;
			} else gauss(tE, tEl, h, 2*m, 2*m);
		}
		Pix* nS = upsample(S, W, H, h, m);
		delete[] S;
		S = nS;
		W *= 2;
		H *= 2;
		if(l < r_size && r[l] > 0)
			jitter(S, W, H, h, m, r[l]);
		if(l > 2) {
			VVP C;
			initCoherence(C, m);
			sprintf(name, "coherence_2_%d.png", l);
			// loadCoherence(C, name);		/* if coherence has been pre_computed then load it */
			// appendCoherence(C, El, m, h);		/* To compute coherence. May take a while */
			// writeCoherence(C, 1, name, m);		/* To save coherence in a file */
			for(int i = 0; i < c; i++)
				correct(S, W, H, h, m, tor ? 2*m : m, C, tor ? tEl : El, kappa);
		}
		sprintf(name, "out_%d.png", l);
		save(S, W, H, El, m, name);
		sprintf(name, "map_%d.png", l);
		saveS(S, W, H, m, name);
	}
	delete tE;
	return S;
}

int main(int argc, char* argv[]) {
	if(argc != 2) {
		cerr << "This is a texture synthetiser." << endl;
		cerr << "You can use it by typing:" << endl;
		cerr << "\t" << argv[0] << " [filename]" << endl;
		cerr << "Where:" << endl;
		cerr << " -filename is the name of the image file used as an example" << endl;
		return 1;
	}
	char* filename = argv[1];
	int m, m2, col;
	uchar* E = stbi_load(filename, &m, &m2, &col, 3);
	if (!E) {
		cerr << "loading failed" << endl;
		return 1;
	}
	bool sq = false;
	if(m != m2) {
		uchar* temp = square(E, m, m2);
		m = min(m, m2);
		free(E);
		E = temp;
		sq = true;
	}
	int W = 3, H = 2;
	VD r = {0.2, 0.3, 0.35, 0.3, 0.2, 0.1, 0.1, 0.1};
	uchar* d = downsample(E, m, 2);
	int ml = m >> 2;
	Pix* S = synthesize(d, ml, r, 3, 0.2, W, H, true);
	uchar* Sh = magnify(ml, E, m, S, W, H);
	stbi_write_png("example.png", m, m, 3, E, 0);
	int mp = pow(2.0, ceil(log2(ml)));
	int Wl = mp*W, Hl = mp*H;
	int Wh = (m*Wl) / ml, Hh = (m*Hl) / ml;
	stbi_write_png("magnific.png", Wh, Hh, 3, Sh, 0);
	delete[] Sh;
	delete[] d;
	if(sq) delete[] E;
	else free(E);
	delete[] S;
	return 0;
}