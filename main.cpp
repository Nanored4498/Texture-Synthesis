#include <iostream>
#include <vector>
#include <queue>
#include <string.h>
#include <sys/stat.h>
#include "filter.hpp"
#include "stb_image.h"
#include "stb_image_write.h"
#include <boost/filesystem.hpp>

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

inline unsigned int int_hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

inline Pix hp(int p) {
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

inline int distNeighb(VI &a, VI &b, int k) {
	int k2 = k*k*3;
	int res = 0, diff;
	for(int i = 0; i < k2; i++) {
		diff = a[i] - b[i];
		res += diff * diff;
	}
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

void appendCoherence(VVP &C, uchar* E, int m, int m2, int h) {
	vector<vector<VI>> Ns(m);
	int kn = min(7, m / h / 2);
	if(kn % 2 == 0) kn ++;
	cout << "KN: " << kn << endl;
	for(int x = 0; x < m; x++) {
		Ns[x] = vector<VI>(m);
		for(int y = 0; y < m; y++)
			getNeighb(x, y, kn, E, m2, Ns[x][y], h);
	}
	double md = pow(0.35*m/kn, 2.0);
	int i = 0;
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
						if(pow(min(dx, m2-dx), 2) + pow(min(dy, m2-dy), 2) < md) {
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
		cerr << "Coherence line: " << y << "/" << m << endl;
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

inline bool is_tor(uchar* im, int m) {
	int d = 0;
	for(int c = 0; c < 3; c++)
	for(int i = 0; i < m; i++)
		d += pow(im[3*i+c] - im[3*(i+m*(m-1))+c], 2) + pow(im[3*i*m+c] - im[3*(m-1+m*i)+c], 2);
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
	int m = min(W, H);
	uchar* res = new uchar[3*m*m];
	for(int c = 0; c < 3; c++)
		for(int x = 0; x < m; x++)
			for(int y = 0; y < m; y++)
				res[3*(x+m*y)+c] = im[3*(x+W*y)+c];
	return res;
}

uchar* downsample(uchar* im, int m, int p) {
	int h = 1 << p;
	int m2 = m >> p;
	uchar* res = new uchar[3*m2*m2];
	double* sum = new double[3*m2*m2];
	for(int i = 0; i < 3*m2*m2; i++) sum[i]=0;
	for(int c = 0; c < 3; c++)
		for(int x = 0; x < m2; x++)
			for(int y = 0; y < m2; y++)
				for(int x2 = h*x; x2 < h*(x+1); x2 ++)
					for(int y2 = h*y; y2 < h*(y+1); y2 ++)
						sum[3*(x+m2*y)+c] += im[3*(x2+m*y2)+c];
	for(int i = 0; i < 3*m2*m2; i++) res[i]=sum[i]/h/h;
	delete[] sum;
	// uchar* g = new uchar[3*m*m];
	// gauss(im, g, h, m, m);
	// for(int c = 0; c < 3; c++)
	// 	for(int x = 0; x < m2; x++)
	// 		for(int y = 0; y < m2; y++)
	// 			res[3*(x+m2*y)+c] = ((int) g[3*h*(x+m*y)+c] + (int) g[3*(h*((x+1)+m*(y+1))-1-m)+c]) >> 1;
	// delete[] g;
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

Pix* synthesize(uchar* E, int m, VD &r, int c, double kappa, int W, int H,
				bool tor=false, const char* file="", bool compute_co=false) {
	// Initialization of images that will be used
	Pix* S = new Pix[W*H];
	int m2 = tor ? 2*m : m;
	uchar* El = new uchar[3*m2*m2];
	uchar *tE = NULL;
	if(tor) tE = torrify(E, m);

	// Folder of coherence
	bool have_folder = strcmp(file, "") != 0;
	char folder[100];
	char name[125];
	if(have_folder) {
		strcpy(folder, file);
		folder[strlen(file)-4] = 0;
	}

	// Some variables
	int L = ceil(log2(m));
	int r_size = r.size();

	// Creation of S_0
	for(int i = 0; i < W*H; i++)
		S[i] = {0, 0};
	if(r_size > 0 && r[0] > 0)
		jitter(S, W, H, m, m, r[0]);
	
	for(int l = 1; l <= L; l++) {
		// footstep
		int h = 1 << (L-l);

		// Computation of gaussian stack
		if(h == 1) {
			delete[] El;
			El = tor ? tE : E;
		} else gauss(tor ? tE : E, El, h*0.7, m2, m2);

		// Upsampling
		Pix* nS = upsample(S, W, H, h, m);
		delete[] S;
		S = nS;
		W *= 2;
		H *= 2;

		// Jitter
		if(l < r_size && r[l] > 0)
			jitter(S, W, H, h, m, r[l]);
		
		// Correction
		if(l > 2) {
			VVP C;
			initCoherence(C, m);
			if(have_folder) {
				sprintf(name, "%s/coherence_%d.png", folder, l);
				struct stat buffer;
				if(stat(name, &buffer) == 0) // Coherence already exists
					loadCoherence(C, name);
				else if(compute_co) { // Coherence does not exists and user want to create it
					boost::filesystem::path dir(folder);
					if(boost::filesystem::create_directory(folder))
						cerr<< "Directory Created: " << folder << endl;
					appendCoherence(C, El, m, m2, h);
					writeCoherence(C, 1, name, m);
				}
			}
			for(int i = 0; i < c; i++)
				correct(S, W, H, h, m, m2, C, El, kappa);
		}

		// Saving images
		sprintf(name, "out_%d.png", l);
		save(S, W, H, El, m2, name);
		sprintf(name, "map_%d.png", l);
		saveS(S, W, H, m, name);
		// sprintf(name, "El_%d.png", l);
		// stbi_write_png(name, m2, m2, 3, El, 0);
	}

	delete[] tE;
	return S;
}

int main(int argc, char* argv[]) {

	// Parameters
	int W = 3, H = 2;
	VD r = {0.2, 0.3, 0.35, 0.3, 0.2, 0.1, 0.1, 0.1};
	bool compute_co = false;

	// Parsing
	if(argc < 2) {
		cerr << "This is a texture synthetiser." << endl;
		cerr << "You can use it by typing:" << endl;
		cerr << "\t" << argv[0] << " [filename] [-c]" << endl;
		cerr << "Where:" << endl;
		cerr << " -filename is the name of the image file used as an example" << endl;
		cerr << " -c if this parameter is given then coherence is computed if it doesn't already exist" << endl;
		return 1;
	}
	for(int i = 2; i < argc; i++) {
		if(strcmp(argv[i], "-c") == 0) compute_co = true;
	}

	// Loading image
	char* filename = argv[1];
	int w, h, col;
	uchar* E = stbi_load(filename, &w, &h, &col, 3);
	if (!E) {
		cerr << "loading failed" << endl;
		return 1;
	}

	// Eventually resize to a square
	bool sq = false;
	int m = w;
	if(w != h) {
		uchar* temp = square(E, w, h);
		m = min(w, h);
		free(E);
		E = temp;
		sq = true;
	}

	// Check if the example is a tore
	bool is_tore = is_tor(E, m);
	if(is_tore)
		cout << "The example has been considered as a tore" << endl;

	// Check if we have to downsize example
	if(m > 256) {
		int ml = m;
		int ds = 0;
		while(ml > 256) ml >>= 1, ds++;
		uchar* d = downsample(E, m, ds);
		Pix* S = synthesize(d, ml, r, 3, 0.2, W, H, !is_tore, filename, compute_co);
		uchar* Sh = magnify(ml, E, m, S, W, H);
		int mp = pow(2.0, ceil(log2(ml)));
		int Wl = mp*W, Hl = mp*H;
		int Wh = (m*Wl) / ml, Hh = (m*Hl) / ml;
		stbi_write_png("magnific.png", Wh, Hh, 3, Sh, 0);
		delete[] S;
		delete[] Sh;
		delete[] d;
	} else {
		Pix* S = synthesize(E, m, r, 3, 0.2, W, H, !is_tore, filename, compute_co);
		int L = 1 << (int) ceil(log2(m));
		save_smooth(S, W*L, H*L, E, m, "magnific.png");
		delete[] S;
	}

	if(sq) delete[] E;
	else free(E);
	return 0;
}