#include "synthetizer.hpp"
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include "filter.hpp"
#include "coherence.hpp"
#include "stb_image.h"
#include "stb_image_write.h"


/****************************/
/*** Phase 1 : Upsampling ***/
/****************************/

Pix* upsample(Pix* S, int W, int H, int h, int md) {
	Pix* res = new Pix[W*H*4];
	int W2 = W*2;
	int h1 = floor(h*0.5);
	int h2 = floor(md - double(h*0.5));
	Pix ze = {h2, h2}, ri = {h1, h2}, up = {h2, h1}, di = {h1, h1};
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			Pix u = S[i+W*j];
			int p = 2*(i + j*W2);
			res[p] = (u + ze) % md;
			res[p+1] = (u + ri) % md;
			res[p+W2] = (u + up) % md;
			res[p+W2+1] = (u + di) % md;
		}
	}
	return res;
}

/************************/
/*** Phase 2 : Jitter ***/
/************************/

inline unsigned int int_hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

inline void hp(int p, double &hpx, double &hpy) {
	unsigned int a = int_hash(p);
	unsigned int b = int_hash(a+p);
	hpx = (double(a % 512) - 256) / 256;
	hpy = (double(b % 512) - 256) / 256;
}

void jitter(Pix* S, int W, int H, int h, int m, double r) {
	double hr = h*r;
	int size = W*H;
	double hpx, hpy;
	for(int i = 0; i < size; i++) {
		hp(i, hpx, hpy);
		S[i] = (S[i] + Pix(m + floor(hr*hpx + 0.5), m + floor(hr*hpy + 0.5))) % m;
	}
}

/****************************/
/*** Phase 3 : Correction ***/
/****************************/

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

/***********************/
/*** Main Algorithm ***/
/**********************/

void init_variables(uchar *Ed, int md, bool tor, const char* file,
					int &m2, uchar*&E2, bool &have_folder, char* folder, int &L) {
	// Initialization of images that will be used
	m2 = tor ? 2*md : md;
	E2 = tor ? torrify(Ed, md) : Ed;
	// Folder of coherence
	have_folder = strcmp(file, "") != 0;
	if(have_folder) {
		strcpy(folder, file);
		folder[strlen(file)-4] = 0;
	}
	// Number of steps
	L = ceil(log2(md));
}

void init_live_WH(int L, int W0, int H0, Pix *S[], int W[], int H[]) {
	W[0] = W0;
	H[0] = H0;
	if(!S[0]) delete[] S[0];
	S[0] = new Pix[W0*H0];
	for(int l = 1; l <= L; l++) {
		W[l] = 2*W[l-1];
		H[l] = 2*H[l-1];
	}
}

void init_live(int W0, int H0, uchar* E2, int m2, int L,
				Pix *S[], int W[], int H[], uchar* El[]) {
	W[0] = W0;
	H[0] = H0;
	S[0] = new Pix[W0*H0];
	for(int l = 1; l <= L; l++) {
		int h = 1 << (L-l);
		if(h == 1)
			El[l] = E2;
		else {
			El[l] = new uchar[3*m2*m2];
			gauss(E2, El[l], h, m2, m2);
		}
		W[l] = 2*W[l-1];
		H[l] = 2*H[l-1];
	}
}

void synthesize_step(int l, Pix *S[], int W[], int H[], uchar *E2, uchar *El[], int md, int m2,
					VD &r, int L, bool have_folder, char* folder, bool compute_co, int c, double kappa,
					bool saveE) {
	if(l == 0) { // Creation of S_0
		for(int i = 0; i < W[0]*H[0]; i++)
			S[0][i] = {0, 0};
		if(r.size() > 0 && r[0] > 0)
			jitter(S[0], W[0], H[0], md, md, r[0]);
		save(S[0], W[0], H[0], E2, m2, "out.png");
	} else {
		char name[100];
		// footstep
		int h = 1 << (L-l);
		// Upsampling
		Pix* nS = upsample(S[l-1], W[l-1], H[l-1], h, md);
		// Jitter
		if(l < (int) r.size() && r[l] > 0)
			jitter(nS, W[l], H[l], h, md, r[l]);
		// Correction
		if(l > 2) {
			VVP C;
			initCoherence(C, md);
			if(have_folder) {
				sprintf(name, "%s/coherence_%d.png", folder, l);
				struct stat buffer;
				if(stat(name, &buffer) == 0) // Coherence already exists
					loadCoherence(C, name);
				else if(compute_co) { // Coherence does not exists and user want to create it
					boost::filesystem::path dir(folder);
					if(boost::filesystem::create_directory(folder))
						std::cerr << "Directory Created: " << folder << std::endl;
					appendCoherence(C, El[l], md, m2, h);
					writeCoherence(C, 1, name, md);
				}
			}
			for(int i = 0; i < c; i++)
				correct(nS, W[l], H[l], h, md, m2, C, El[l], kappa);
		}
		if(saveE)
			save(nS, W[l], H[l], El[l], m2, "out.png");
		else
			saveS(nS, W[l], H[l], md, "out.png");
		if(!S[l]) delete[] S[l];
		S[l] = nS;
	}
}

Pix* synthesize(uchar* E, int m, VD &r, int c, double kappa, int W, int H,
				bool tor, const char* file, bool compute_co) {
	// Initialization of variables
	Pix* S = new Pix[W*H];
	int m2;
	uchar *E2;
	bool have_folder;
	char folder[100], name[125];
	int L;
	init_variables(E, m, tor, file,
					m2, E2, have_folder, folder, L);
	uchar* El = new uchar[3*m2*m2];

	// Creation of S_0
	for(int i = 0; i < W*H; i++)
		S[i] = {0, 0};
	if(r.size() > 0 && r[0] > 0)
		jitter(S, W, H, m, m, r[0]);
	
	for(int l = 1; l <= L; l++) {
		// footstep
		int h = 1 << (L-l);

		// Computation of gaussian stack
		if(h == 1) {
			delete[] El;
			El = E2;
		} else gauss(E2, El, h, m2, m2);

		// Upsampling
		Pix* nS = upsample(S, W, H, h, m);
		delete[] S;
		S = nS;
		W *= 2;
		H *= 2;

		// Jitter
		if(l < (int) r.size() && r[l] > 0)
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
						std::cerr << "Directory Created: " << folder << std::endl;
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

	if(tor) delete[] E2;
	return S;
}

/**********************/
/*** Magnification ***/
/*********************/

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
	return res;
}

uchar* magnify(int ml, uchar* Eh, int mh, Pix* S, int W, int H, int &Wh, int &Hh) {
	int mp = pow(2.0, ceil(log2(ml)));
	int Wl = mp*W, Hl = mp*H;
	Wh = (mh*Wl) / ml, Hh = (mh*Hl) / ml;
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
							colors[i][j] = Eh[3 * (std::max(0, u.x) + mh*std::max(0, u.y)) + c];
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

/***********************************/
/*** Load Image and pretreatment ***/
/***********************************/

int load_image(const char* filename, double to_tor, uchar *&E, int &m, double &is_tore, double &new_E, uchar *&Ed, int &md) {
	int w, h, col;
	E = stbi_load(filename, &w, &h, &col, 3);
	if (!E) {
		std::cerr << "loading failed" << std::endl;
		return 1;
	}

	// Eventually resize to a square
	new_E = false;
	m = std::min(w, h);
	if(w != h) {
		uchar* temp = square(E, w, h);
		free(E);
		E = temp;
		new_E = true;
	}

	// If the user want to torify image then torify the image
	if(to_tor) {
		uchar* tE = torrify(E, m);
		if(new_E) delete[] E;
		else free(E);
		E = tE;
		m *= 2;
		new_E = true;
	}

	// Check if the example is a tore
	is_tore = is_tor(E, m);
	if(is_tore)
		std::cerr << "The example has been considered as a tore" << std::endl;

	// Downsizing
	md = m;
	if(m > 256) {
		int ds = 0;
		while(md > 256) md >>= 1, ds++;
		Ed = downsample(E, m, ds);
	} else Ed = E;
	return 0;
}